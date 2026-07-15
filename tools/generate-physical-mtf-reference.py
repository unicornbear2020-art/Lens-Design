#!/usr/bin/env python3
"""Generate pinned Optiland physical-MTF coverage fixtures.

Requires optiland==0.5.8 and numpy. The generated JavaScript is runtime-only
data; Lens Design never imports Optiland or performs network access.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
from pathlib import Path

import numpy as np
from optiland.mtf.fft import FFTMTF
from optiland.psf.fft import FFTPSF
from optiland.samples.objectives import DoubleGauss, ReverseTelephoto


FREQUENCIES_LP_MM = (0, 10, 30, 40, 50)
DOUBLE_GAUSS_FIELDS = (
    {"key": "center", "name": "On-axis", "coordinates": (0, 0), "angleDeg": 0},
    {"key": "mid", "name": "Mid field", "coordinates": (0, 10 / 14), "angleDeg": 10},
    {"key": "corner", "name": "Full field", "coordinates": (0, 1), "angleDeg": 14},
)
REVERSE_TELEPHOTO_FIELDS = (
    {"key": "center", "name": "On-axis", "coordinates": (0, 0), "angleDeg": 0},
    {"key": "mid", "name": "Mid field", "coordinates": (0, 21 / 30), "angleDeg": 21},
    {"key": "corner", "name": "Full field", "coordinates": (0, 1), "angleDeg": 30},
)
WAVELENGTHS = (
    {"key": "C", "wavelengthUm": 0.6563},
    {"key": "d", "wavelengthUm": 0.5876},
    {"key": "F", "wavelengthUm": 0.4861},
)
REQUESTED_RAYS = 256
PUPIL_RAYS = 90
FFT_GRID_SIZE = 512


def encode_pupil(pupil: np.ndarray) -> dict:
    amplitude = np.abs(pupil)
    phase = np.angle(pupil)
    mask = amplitude > 0
    amplitude_quantized = np.zeros(amplitude.shape, dtype=np.uint8)
    amplitude_quantized[mask] = np.clip(
        np.rint((amplitude[mask] - 1) * 1_000_000) + 128,
        1,
        255,
    ).astype(np.uint8)
    phase_quantized = np.zeros(phase.shape, dtype="<i2")
    phase_quantized[mask] = np.rint(phase[mask] / np.pi * 32767).astype("<i2")
    payload = amplitude_quantized.tobytes() + phase_quantized.tobytes()
    return {
        "size": PUPIL_RAYS,
        "amplitude": "uint8, 1 ppm steps around unity; 0 is outside pupil",
        "phase": "little-endian int16 mapped from [-pi,+pi] to [-32767,+32767]",
        "payloadSha256": hashlib.sha256(payload).hexdigest(),
        "nonzeroSamples": int(mask.sum()),
        "amplitudeBase64": base64.b64encode(amplitude_quantized.tobytes()).decode(),
        "phaseBase64": base64.b64encode(phase_quantized.tobytes()).decode(),
    }


def build_case(optic, lens_slug: str, field: dict, wavelength: dict) -> dict:
    coordinates = field["coordinates"]
    wavelength_um = wavelength["wavelengthUm"]
    mtf = FFTMTF(
        optic,
        fields=[coordinates],
        wavelength=wavelength_um,
        num_rays=REQUESTED_RAYS,
        max_freq=100,
        strategy="chief_ray",
        remove_tilt=False,
    )
    psf = FFTPSF(
        optic,
        coordinates,
        wavelength_um,
        num_rays=PUPIL_RAYS,
        grid_size=FFT_GRID_SIZE,
        strategy="chief_ray",
        remove_tilt=False,
    )

    frequency = np.asarray(mtf.freq)
    tangential = np.asarray(mtf.mtf[0][0])
    sagittal = np.asarray(mtf.mtf[0][1])
    samples = [
        {
            "frequencyLpMm": value,
            "tangential": float(np.interp(value, frequency, tangential)),
            "sagittal": float(np.interp(value, frequency, sagittal)),
        }
        for value in FREQUENCIES_LP_MM
    ]

    return {
        "id": f"optiland-v0.5.8-{lens_slug}-{field['key']}-{wavelength['key']}",
        "fieldKey": field["key"],
        "fieldName": field["name"],
        "field": list(coordinates),
        "fieldAngleDeg": field["angleDeg"],
        "wavelengthKey": wavelength["key"],
        "wavelengthUm": wavelength_um,
        "workingFNumber": float(mtf.FNO),
        "frequencyStepLpMm": float(frequency[1] - frequency[0]),
        "effectivePupilRays": int(mtf.num_rays),
        "pupilEncoding": encode_pupil(np.asarray(psf.pupils[0])),
        "samples": samples,
        "tolerance": 0.005,
    }


def build_fixture(
    optic,
    *,
    fixture_id: str,
    lens_slug: str,
    lens_class: str,
    fields: tuple,
    source_url: str,
    scope: str,
) -> dict:
    cases = [
        build_case(optic, lens_slug, field, wavelength)
        for field in fields
        for wavelength in WAVELENGTHS
    ]
    return {
        "id": fixture_id,
        "schemaVersion": 2,
        "solver": "Optiland FFTMTF",
        "solverVersion": "0.5.8",
        "sourceUrl": source_url,
        "sourceBlobSha": "69044dae6866b187cc627b74345804af393ae62f",
        "lensClass": lens_class,
        "strategy": "chief_ray",
        "removeTilt": False,
        "requestedNumRays": REQUESTED_RAYS,
        "fftGridSize": FFT_GRID_SIZE,
        "fieldType": "angle",
        "fields": [
            {key: value for key, value in field.items() if key != "coordinates"}
            | {"coordinates": list(field["coordinates"])}
            for field in fields
        ],
        "wavelengths": list(WAVELENGTHS),
        "frequenciesLpMm": list(FREQUENCIES_LP_MM),
        "cases": cases,
        "tolerance": 0.005,
        "coverage": {
            "fieldCount": len(fields),
            "wavelengthCount": len(WAVELENGTHS),
            "caseCount": len(cases),
            "comparisonCount": len(cases) * len(FREQUENCIES_LP_MM) * 2,
        },
        "scope": scope,
    }


def build_fixtures() -> tuple[dict, dict]:
    double_gauss = build_fixture(
        DoubleGauss(),
        fixture_id="optiland-v0.5.8-double-gauss-field-spectral-matrix",
        lens_slug="double-gauss",
        lens_class="DoubleGauss",
        fields=DOUBLE_GAUSS_FIELDS,
        source_url="https://github.com/optiland/optiland/blob/v0.5.8/optiland/samples/objectives.py#L75-L114",
        scope="on-axis and off-axis 0/10/14 degree fields, monochromatic C/d/F, complex exit-pupil to MTF parity",
    )
    reverse_telephoto = build_fixture(
        ReverseTelephoto(),
        fixture_id="optiland-v0.5.8-reverse-telephoto-field-spectral-matrix",
        lens_slug="reverse-telephoto",
        lens_class="ReverseTelephoto",
        fields=REVERSE_TELEPHOTO_FIELDS,
        source_url="https://github.com/optiland/optiland/blob/v0.5.8/optiland/samples/objectives.py#L117-L175",
        scope="independent retrofocus validation at on-axis and off-axis 0/21/30 degree fields, monochromatic C/d/F, complex exit-pupil to MTF parity",
    )
    return double_gauss, reverse_telephoto


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    double_gauss, reverse_telephoto = build_fixtures()
    primary_fixture = json.dumps(double_gauss, separators=(",", ":"))
    secondary_fixture = json.dumps(reverse_telephoto, separators=(",", ":"))
    args.output.write_text(
        "globalThis.PHYSICAL_MTF_EXTERNAL_REFERENCE_FIXTURE = Object.freeze("
        f"{primary_fixture});\n"
        "globalThis.PHYSICAL_MTF_SECONDARY_EXTERNAL_REFERENCE_FIXTURE = Object.freeze("
        f"{secondary_fixture});\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()

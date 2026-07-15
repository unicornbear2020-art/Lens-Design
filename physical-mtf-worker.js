/* eslint-env worker */

const PHYSICAL_MTF_WORKER_VERSION = "20260715-ray-sampling-performance-1";
const PHYSICAL_LENS_PUPIL_RADIUS_FRACTION = 0.44;
const clamp = (value, min, max) => Math.min(max, Math.max(min, Number.isFinite(Number(value)) ? Number(value) : min));
const isPowerOfTwo = (value) => Number.isInteger(value) && value > 0 && (value & (value - 1)) === 0;

const fftRadix2Complex = (inputReal, inputImag = [], inverse = false) => {
  const n = inputReal.length;
  if (n === 0 || (n & (n - 1)) !== 0) return { status: "invalid", warning: "FFT length must be a power of two." };
  const real = inputReal.slice();
  const imag = Array.from({ length: n }, (_, index) => inputImag[index] || 0);
  let j = 0;
  for (let i = 1; i < n; i += 1) {
    let bit = n >> 1;
    while (j & bit) {
      j ^= bit;
      bit >>= 1;
    }
    j ^= bit;
    if (i < j) {
      [real[i], real[j]] = [real[j], real[i]];
      [imag[i], imag[j]] = [imag[j], imag[i]];
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const angle = (inverse ? 2 : -2) * Math.PI / len;
    const wLenReal = Math.cos(angle);
    const wLenImag = Math.sin(angle);
    for (let i = 0; i < n; i += len) {
      let wReal = 1;
      let wImag = 0;
      for (let k = 0; k < len / 2; k += 1) {
        const even = i + k;
        const odd = even + len / 2;
        const oddReal = real[odd] * wReal - imag[odd] * wImag;
        const oddImag = real[odd] * wImag + imag[odd] * wReal;
        real[odd] = real[even] - oddReal;
        imag[odd] = imag[even] - oddImag;
        real[even] += oddReal;
        imag[even] += oddImag;
        const nextReal = wReal * wLenReal - wImag * wLenImag;
        wImag = wReal * wLenImag + wImag * wLenReal;
        wReal = nextReal;
      }
    }
  }
  if (inverse) {
    for (let i = 0; i < n; i += 1) {
      real[i] /= n;
      imag[i] /= n;
    }
  }
  return { status: "valid", real, imag };
};

const fft2DRadix2Complex = (inputReal, inputImag, width, height, inverse = false) => {
  const real = inputReal.slice();
  const imag = inputImag.slice();
  for (let y = 0; y < height; y += 1) {
    const start = y * width;
    const row = fftRadix2Complex(real.slice(start, start + width), imag.slice(start, start + width), inverse);
    if (row.status !== "valid") return row;
    for (let x = 0; x < width; x += 1) {
      real[start + x] = row.real[x];
      imag[start + x] = row.imag[x];
    }
  }
  for (let x = 0; x < width; x += 1) {
    const columnReal = [];
    const columnImag = [];
    for (let y = 0; y < height; y += 1) {
      const index = y * width + x;
      columnReal.push(real[index]);
      columnImag.push(imag[index]);
    }
    const column = fftRadix2Complex(columnReal, columnImag, inverse);
    if (column.status !== "valid") return column;
    for (let y = 0; y < height; y += 1) {
      const index = y * width + x;
      real[index] = column.real[y];
      imag[index] = column.imag[y];
    }
  }
  return { status: "valid", real, imag };
};

const buildIdealCircularPupil = (gridSize) => {
  const size = [128, 256].includes(Number(gridSize)) ? Number(gridSize) : 128;
  const radius = size * 0.16;
  const center = (size - 1) / 2;
  const real = new Array(size * size).fill(0);
  const imag = new Array(size * size).fill(0);
  for (let y = 0; y < size; y += 1) {
    for (let x = 0; x < size; x += 1) {
      if (Math.hypot(x - center, y - center) <= radius) real[y * size + x] = 1;
    }
  }
  return { gridSize: size, radiusPx: radius, real, imag };
};

const analyticMtf = (normalizedFrequency) => {
  const nu = clamp(normalizedFrequency, 0, 1);
  if (nu >= 1) return 0;
  return (2 / Math.PI) * (Math.acos(nu) - nu * Math.sqrt(Math.max(0, 1 - nu ** 2)));
};

const samplePupilGridBilinear = (pupil, x, y) => {
  const size = pupil.gridSize;
  if (!Number.isFinite(x) || !Number.isFinite(y) || x < 0 || y < 0 || x >= size - 1 || y >= size - 1) return 0;
  const x0 = Math.floor(x);
  const y0 = Math.floor(y);
  const dx = x - x0;
  const dy = y - y0;
  const i00 = y0 * size + x0;
  const i10 = i00 + 1;
  const i01 = i00 + size;
  const i11 = i01 + 1;
  const top = pupil.real[i00] * (1 - dx) + pupil.real[i10] * dx;
  const bottom = pupil.real[i01] * (1 - dx) + pupil.real[i11] * dx;
  return top * (1 - dy) + bottom * dy;
};

const sampleComplexPupilGridBilinear = (pupil, x, y) => {
  const size = pupil.gridSize;
  if (!Number.isFinite(x) || !Number.isFinite(y) || x < 0 || y < 0 || x >= size - 1 || y >= size - 1) {
    return { real: 0, imag: 0 };
  }
  const x0 = Math.floor(x);
  const y0 = Math.floor(y);
  const tx = x - x0;
  const ty = y - y0;
  const i00 = y0 * size + x0;
  const i10 = i00 + 1;
  const i01 = i00 + size;
  const i11 = i01 + 1;
  const interpolate = (values) => {
    const top = values[i00] * (1 - tx) + values[i10] * tx;
    const bottom = values[i01] * (1 - tx) + values[i11] * tx;
    return top * (1 - ty) + bottom * ty;
  };
  return { real: interpolate(pupil.real), imag: interpolate(pupil.imag) };
};

const complexPupilAutocorrelationOtf = (pupil, shiftX, shiftY, pupilEnergy) => {
  if (!(pupilEnergy > 0)) return { real: NaN, imag: NaN, magnitude: NaN };
  let correlationReal = 0;
  let correlationImag = 0;
  for (let y = 0; y < pupil.gridSize; y += 1) {
    for (let x = 0; x < pupil.gridSize; x += 1) {
      const index = y * pupil.gridSize + x;
      const real = pupil.real[index] || 0;
      const imag = pupil.imag[index] || 0;
      if (Math.abs(real) + Math.abs(imag) < 1e-15) continue;
      const shifted = sampleComplexPupilGridBilinear(pupil, x + shiftX, y + shiftY);
      correlationReal += real * shifted.real + imag * shifted.imag;
      correlationImag += imag * shifted.real - real * shifted.imag;
    }
  }
  const real = correlationReal / pupilEnergy;
  const imag = correlationImag / pupilEnergy;
  return { real, imag, magnitude: clamp(Math.hypot(real, imag), 0, 1) };
};

const complexPupilAutocorrelationMtf = (pupil, shiftX, shiftY, pupilEnergy) => (
  complexPupilAutocorrelationOtf(pupil, shiftX, shiftY, pupilEnergy).magnitude
);

const sampleComplexPupilOtfAxis = (pupil, offset, dx, dy, pupilEnergy) => {
  const values = [
    complexPupilAutocorrelationMtf(pupil, offset * dx, offset * dy, pupilEnergy),
    complexPupilAutocorrelationMtf(pupil, -offset * dx, -offset * dy, pupilEnergy)
  ].filter(Number.isFinite);
  return values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : NaN;
};

const sampleComplexPupilComplexOtfAxis = (pupil, offset, dx, dy, pupilEnergy) => (
  complexPupilAutocorrelationOtf(pupil, offset * dx, offset * dy, pupilEnergy)
);

const autocorrelationMtf = (pupil, normalizedFrequency) => {
  const nu = clamp(normalizedFrequency, 0, 1);
  if (nu <= 0) return 1;
  if (nu >= 1) return 0;
  const shift = nu * 2 * pupil.radiusPx;
  let total = 0;
  let horizontalOverlap = 0;
  let verticalOverlap = 0;
  for (let y = 0; y < pupil.gridSize; y += 1) {
    for (let x = 0; x < pupil.gridSize; x += 1) {
      const value = pupil.real[y * pupil.gridSize + x];
      if (!value) continue;
      total += value;
      horizontalOverlap += value * samplePupilGridBilinear(pupil, x + shift, y);
      verticalOverlap += value * samplePupilGridBilinear(pupil, x, y + shift);
    }
  }
  return clamp(((horizontalOverlap + verticalOverlap) / 2) / Math.max(total, 1), 0, 1);
};

const sampleOtfMagnitudeUnshifted = (otf, x, y) => {
  const size = otf.gridSize;
  const wrap = (value) => {
    const wrapped = value % size;
    return wrapped < 0 ? wrapped + size : wrapped;
  };
  const sx = wrap(x);
  const sy = wrap(y);
  const x0 = Math.floor(sx);
  const y0 = Math.floor(sy);
  const x1 = (x0 + 1) % size;
  const y1 = (y0 + 1) % size;
  const dx = sx - x0;
  const dy = sy - y0;
  const mag = (px, py) => {
    const index = py * size + px;
    return Math.hypot(otf.real[index], otf.imag[index]);
  };
  const top = mag(x0, y0) * (1 - dx) + mag(x1, y0) * dx;
  const bottom = mag(x0, y1) * (1 - dx) + mag(x1, y1) * dx;
  return top * (1 - dy) + bottom * dy;
};

const sampleOtfRadialMtf = (otf, pupil, normalizedFrequency) => {
  const nu = clamp(normalizedFrequency, 0, 1);
  if (nu <= 0) return 1;
  if (nu >= 1) return 0;
  const radius = nu * 2 * pupil.radiusPx;
  const directions = [[1, 0], [-1, 0], [0, 1], [0, -1], [Math.SQRT1_2, Math.SQRT1_2], [-Math.SQRT1_2, Math.SQRT1_2], [Math.SQRT1_2, -Math.SQRT1_2], [-Math.SQRT1_2, -Math.SQRT1_2]];
  const values = directions.map(([dx, dy]) => sampleOtfMagnitudeUnshifted(otf, radius * dx, radius * dy)).filter(Number.isFinite);
  return values.length ? clamp(values.reduce((sum, value) => sum + value, 0) / values.length, 0, 1) : NaN;
};

const sampleOtfAxisMtf = (otf, offset, dx, dy) => {
  const values = [
    sampleOtfMagnitudeUnshifted(otf, offset * dx, offset * dy),
    sampleOtfMagnitudeUnshifted(otf, -offset * dx, -offset * dy)
  ].filter(Number.isFinite);
  return values.length ? clamp(values.reduce((sum, value) => sum + value, 0) / values.length, 0, 1) : NaN;
};

const calculateCore = (gridSize) => {
  const started = performance.now();
  const pupil = buildIdealCircularPupil(gridSize);
  const field = fft2DRadix2Complex(pupil.real, pupil.imag, pupil.gridSize, pupil.gridSize, false);
  if (field.status !== "valid") return field;
  const psf = field.real.map((real, index) => real ** 2 + field.imag[index] ** 2);
  const psfEnergy = psf.reduce((sum, value) => sum + value, 0) || 1;
  const normalizedPsf = psf.map((value) => value / psfEnergy);
  const otf = fft2DRadix2Complex(normalizedPsf, new Array(normalizedPsf.length).fill(0), pupil.gridSize, pupil.gridSize, false);
  if (otf.status !== "valid") return otf;
  const dcReal = otf.real[0] || 0;
  const dcImag = otf.imag[0] || 0;
  const dc = Math.hypot(dcReal, dcImag) || 1;
  const normalizedOtf = { gridSize: pupil.gridSize, real: otf.real.map((value) => value / dc), imag: otf.imag.map((value) => value / dc) };
  const samples = Array.from({ length: 101 }, (_, index) => {
    const normalizedFrequency = index / 100;
    const fftValue = sampleOtfRadialMtf(normalizedOtf, pupil, normalizedFrequency);
    const analyticValue = analyticMtf(normalizedFrequency);
    const autocorrelationValue = autocorrelationMtf(pupil, normalizedFrequency);
    return { normalizedFrequency, fftValue, analyticValue, autocorrelationValue, absoluteError: Math.abs(fftValue - analyticValue) };
  });
  const maxError = Math.max(...samples.map((sample) => sample.absoluteError));
  const rmsError = Math.sqrt(samples.reduce((sum, sample) => sum + sample.absoluteError ** 2, 0) / samples.length);
  const tolerance = pupil.gridSize >= 256 ? { maxError: 0.01, rmsError: 0.005 } : { maxError: 0.02, rmsError: 0.01 };
  const autocorrelationErrors = samples.map((sample) => Math.abs(sample.autocorrelationValue - sample.analyticValue));
  const fftAutocorrelationErrors = samples.map((sample) => Math.abs(sample.fftValue - sample.autocorrelationValue));
  const maxAutocorrelationError = Math.max(...autocorrelationErrors);
  const rmsAutocorrelationError = Math.sqrt(autocorrelationErrors.reduce((sum, error) => sum + error ** 2, 0) / autocorrelationErrors.length);
  const mtfZero = samples[0]?.fftValue ?? NaN;
  const mtfCutoff = samples[samples.length - 1]?.fftValue ?? NaN;
  const sanityPassed = Math.abs(mtfZero - 1) < 0.001 && mtfCutoff <= 0.05 && samples.every((sample) => sample.fftValue >= -1e-9);
  const matchPassed = maxError <= tolerance.maxError && rmsError <= tolerance.rmsError;
  return {
    status: matchPassed && sanityPassed ? "fft-match-passed" : "prototype-failed",
    validationLabel: matchPassed && sanityPassed ? "FFT pipeline match passed" : "Prototype — FFT validation failed",
    gridSize: pupil.gridSize,
    samples,
    fftDerivedMtf: samples.map(({ normalizedFrequency, fftValue }) => ({ normalizedFrequency, value: fftValue })),
    analyticMtf: samples.map(({ normalizedFrequency, analyticValue }) => ({ normalizedFrequency, value: analyticValue })),
    autocorrelationReferenceMtf: samples.map(({ normalizedFrequency, autocorrelationValue }) => ({ normalizedFrequency, value: autocorrelationValue })),
    maxError,
    rmsError,
    maxAutocorrelationError,
    rmsAutocorrelationError,
    maxFftAutocorrelationDifference: Math.max(...fftAutocorrelationErrors),
    elapsedMs: performance.now() - started,
    tolerance,
    validation: {
      autocorrelationReference: { maxError: maxAutocorrelationError, rmsError: rmsAutocorrelationError, passed: maxAutocorrelationError <= tolerance.maxError && rmsAutocorrelationError <= tolerance.rmsError },
      fftPipeline: { maxError, rmsError, passed: matchPassed && sanityPassed },
      convergence: null
    },
    pipeline: {
      psfEnergy,
      normalizedPsfEnergy: normalizedPsf.reduce((sum, value) => sum + value, 0),
      otfDc: dc,
      otfDcReal: dcReal,
      otfDcImag: dcImag,
      mtfZero,
      mtfCutoff,
      note: "Worker FFT-derived curve is sampled from OTF generated by complex pupil → FFT PSF → FFT OTF."
    }
  };
};

const calculateReferencePupilParity = (pupil, cutoffFrequencyLpMm) => {
  const amplitude = pupil.real.map((real, index) => Math.hypot(real, pupil.imag[index] || 0));
  const referencePupil = {
    gridSize: pupil.gridSize,
    radiusPx: pupil.radiusPx,
    real: amplitude,
    imag: new Array(amplitude.length).fill(0)
  };
  const center = (pupil.gridSize - 1) / 2;
  let circularPixelCount = 0;
  let transmittedPixelCount = 0;
  let amplitudeSum = 0;
  amplitude.forEach((value, index) => {
    const x = index % pupil.gridSize;
    const y = Math.floor(index / pupil.gridSize);
    if (Math.hypot(x - center, y - center) > pupil.radiusPx) return;
    circularPixelCount += 1;
    if (value > 0.001) transmittedPixelCount += 1;
    amplitudeSum += value;
  });
  const pupilFillFraction = circularPixelCount > 0 ? transmittedPixelCount / circularPixelCount : 0;
  const meanAmplitude = circularPixelCount > 0 ? amplitudeSum / circularPixelCount : 0;
  const applicable = pupilFillFraction >= 0.98 && meanAmplitude >= 0.95 && cutoffFrequencyLpMm > 0;
  const energy = amplitude.reduce((sum, value) => sum + value ** 2, 0);
  const comparisons = applicable
    ? [10, 30, 40, 50].flatMap((frequencyLpMm) => {
      const normalizedFrequency = clamp(frequencyLpMm / cutoffFrequencyLpMm, 0, 1);
      const offset = normalizedFrequency * 2 * pupil.radiusPx;
      const analyticValue = analyticMtf(normalizedFrequency);
      return [
        ["sagittal", 0, 1],
        ["tangential", 1, 0]
      ].map(([axis, dx, dy]) => {
        const referenceValue = sampleComplexPupilOtfAxis(referencePupil, offset, dx, dy, energy);
        return {
          frequencyLpMm,
          normalizedFrequency,
          axis,
          referenceValue,
          analyticValue,
          absoluteDifference: Math.abs(referenceValue - analyticValue)
        };
      });
    })
    : [];
  const tolerance = pupil.gridSize >= 256 ? 0.02 : 0.035;
  const maximumDifference = comparisons.length
    ? Math.max(...comparisons.map((item) => item.absoluteDifference))
    : NaN;
  const passed = applicable
    && comparisons.length === 8
    && comparisons.every((item) => Number.isFinite(item.absoluteDifference) && item.absoluteDifference <= tolerance);
  return {
    status: !applicable ? "not-applicable" : passed ? "passed" : "failed",
    method: "zero-phase reference pupil compared with analytic circular-aperture MTF",
    applicable,
    passed,
    tolerance,
    maximumDifference,
    pupilFillFraction,
    meanAmplitude,
    comparisons
  };
};

const calculateLensComplexPupilMtf = (payload = {}) => {
  const started = performance.now();
  const gridSize = Number(payload.gridSize);
  const expectedLength = gridSize * gridSize;
  const real = payload.realBuffer instanceof ArrayBuffer
    ? Array.from(new Float64Array(payload.realBuffer))
    : Array.from(payload.real || []);
  const imag = payload.imagBuffer instanceof ArrayBuffer
    ? Array.from(new Float64Array(payload.imagBuffer))
    : Array.from(payload.imag || []);
  if (!isPowerOfTwo(gridSize) || real.length !== expectedLength || imag.length !== expectedLength) {
    return {
      status: "invalid",
      reason: "Complex pupil arrays must match a power-of-two square grid."
    };
  }
  const radiusPx = Number(payload.radiusPx) > 0
    ? Number(payload.radiusPx)
    : gridSize * PHYSICAL_LENS_PUPIL_RADIUS_FRACTION;
  const wavelengthNm = Number(payload.wavelengthNm);
  const workingFNumber = Number(payload.workingFNumber);
  const wavelengthMm = wavelengthNm / 1000000;
  const cutoffFrequencyLpMm = wavelengthMm > 0 && workingFNumber > 0
    ? 1 / (wavelengthMm * workingFNumber)
    : NaN;
  const pupil = { gridSize, radiusPx, real, imag };
  const pupilEnergy = real.reduce((sum, value, index) => sum + value ** 2 + (imag[index] || 0) ** 2, 0);
  const field = fft2DRadix2Complex(real, imag, gridSize, gridSize, false);
  if (field.status !== "valid") return field;
  const psf = field.real.map((value, index) => value ** 2 + field.imag[index] ** 2);
  const psfEnergy = psf.reduce((sum, value) => sum + value, 0);
  if (!(psfEnergy > 0)) return { status: "invalid", reason: "Complex pupil has no transmitted energy." };
  const normalizedPsf = psf.map((value) => value / psfEnergy);
  const otfResult = fft2DRadix2Complex(normalizedPsf, new Array(expectedLength).fill(0), gridSize, gridSize, false);
  if (otfResult.status !== "valid") return otfResult;
  const otfDcReal = otfResult.real[0] || 0;
  const otfDcImag = otfResult.imag[0] || 0;
  const otfDc = Math.hypot(otfDcReal, otfDcImag) || 1;
  const otf = {
    gridSize,
    real: otfResult.real.map((value) => value / otfDc),
    imag: otfResult.imag.map((value) => value / otfDc)
  };
  const nativeFrequencyStepLpMm = Number(payload.frequencyStepLpMm);
  const usesNativePupilFrequencyGrid = nativeFrequencyStepLpMm > 0;
  const sampleCount = usesNativePupilFrequencyGrid
    ? Math.max(21, Math.floor(2 * radiusPx + 0.000001) + 1)
    : Math.max(21, Math.round(Number(payload.sampleCount) || 101));
  const makeCurve = (dx, dy) => Array.from({ length: sampleCount }, (_, index) => {
    const pupilOffset = usesNativePupilFrequencyGrid ? index : index / (sampleCount - 1) * 2 * radiusPx;
    const normalizedFrequency = pupilOffset / (2 * radiusPx);
    const complexOtf = normalizedFrequency >= 1
      ? { real: 0, imag: 0, magnitude: 0 }
      : sampleComplexPupilComplexOtfAxis(pupil, pupilOffset, dx, dy, pupilEnergy);
    return {
      normalizedFrequency,
      lpMm: usesNativePupilFrequencyGrid
        ? index * nativeFrequencyStepLpMm
        : Number.isFinite(cutoffFrequencyLpMm) ? normalizedFrequency * cutoffFrequencyLpMm : NaN,
      value: complexOtf.magnitude,
      otfReal: complexOtf.real,
      otfImag: complexOtf.imag,
      fftGridValue: normalizedFrequency >= 1
        ? 0
        : sampleOtfAxisMtf(otf, pupilOffset, dx, dy)
    };
  });
  // Standard fields vary along image-space Y: tangential is the meridional Y cut,
  // while sagittal is the orthogonal X cut. Keep this identical to app.js fallback.
  const tangential = makeCurve(0, 1);
  const sagittal = makeCurve(1, 0);
  const finiteValues = [...tangential, ...sagittal].map((point) => point.value).filter(Number.isFinite);
  const zeroSanityPassed = Math.abs((tangential[0]?.value ?? NaN) - 1) < 0.000001
    && Math.abs((sagittal[0]?.value ?? NaN) - 1) < 0.000001;
  const rangeSanityPassed = finiteValues.length === sampleCount * 2
    && finiteValues.every((value) => value >= -1e-9 && value <= 1.000001);
  const referenceParity = calculateReferencePupilParity(pupil, cutoffFrequencyLpMm);
  return {
    status: zeroSanityPassed && rangeSanityPassed ? "opd-connected-worker-prototype" : "prototype-failed",
    phase: "Phase 2 Step 4",
    gridSize,
    radiusPx,
    wavelengthNm,
    workingFNumber,
    cutoffFrequencyLpMm,
    mtf: { tangential, sagittal },
    referenceParity,
    elapsedMs: performance.now() - started,
    validation: {
      zeroFrequency: { passed: zeroSanityPassed },
      range: { passed: rangeSanityPassed },
      convergence: null
    },
    diagnostics: {
      normalizedPsfEnergy: normalizedPsf.reduce((sum, value) => sum + value, 0),
      otfDc,
      otfDcReal,
      otfDcImag,
      mtfZeroTangential: tangential[0]?.value,
      mtfZeroSagittal: sagittal[0]?.value,
      minimumMtf: finiteValues.length ? Math.min(...finiteValues) : NaN,
      maximumMtf: finiteValues.length ? Math.max(...finiteValues) : NaN,
      cutSamplingMethod: "subpixel complex-pupil autocorrelation (FFT-equivalent OTF)"
    }
  };
};

self.onmessage = (event) => {
  const payload = event.data || {};
  const { requestId, gridSize, task } = payload;
  try {
    const result = task === "lens-complex-pupil"
      ? calculateLensComplexPupilMtf(payload)
      : calculateCore(gridSize);
    self.postMessage({
      requestId,
      status: result.status === "invalid" || result.status === "prototype-failed" ? "error" : "complete",
      workerVersion: PHYSICAL_MTF_WORKER_VERSION,
      task: task || "ideal-circular-pupil",
      result
    });
  } catch (error) {
    self.postMessage({
      requestId,
      status: "error",
      workerVersion: PHYSICAL_MTF_WORKER_VERSION,
      task: task || "ideal-circular-pupil",
      result: { status: "invalid", warning: error?.message || "Worker calculation failed." }
    });
  }
};

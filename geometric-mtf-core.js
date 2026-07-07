/* eslint-env worker */

(() => {
  const GEOMETRIC_MTF_CORE_VERSION = "20260707-ipad-system-result-1";
  const GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION = "geometric-lsf-contract-20260630-1";
  const QUALITY_PROFILES = {
    interactive: { baseGrid: 16, label: "Interactive" },
    high: { baseGrid: 32, label: "High" },
    reference: { baseGrid: 64, label: "Reference" }
  };
  const MAIN_PANEL_QUALITY_PROFILES = {
    interactive: { rings: 6, label: "Interactive", approximateRays: 129 },
    high: { rings: 15, label: "High", approximateRays: 723 },
    reference: { rings: 20, label: "Reference", approximateRays: 1263 }
  };
  const MTF_READOUT_FREQUENCIES = [10, 20, 40, 80];

  const toNumber = (value) => {
    const number = Number(value);
    return Number.isFinite(number) ? number : NaN;
  };
  const clamp = (value, min, max) => Math.min(max, Math.max(min, value));
  const normalizeVector3 = (x, y, z) => {
    const length = Math.hypot(x, y, z);
    return length > 1e-12 ? { x: x / length, y: y / length, z: z / length } : null;
  };
  const dot3 = (a, b) => a.x * b.x + a.y * b.y + a.z * b.z;
  const isPlanoRadius = (radius) => !Number.isFinite(toNumber(radius)) || Math.abs(toNumber(radius)) < 1e-9;
  const nextPowerOfTwo = (value) => {
    let power = 1;
    while (power < value) power *= 2;
    return power;
  };
  const sinc = (value) => Math.abs(value) < 1e-12 ? 1 : Math.sin(Math.PI * value) / (Math.PI * value);
  const asphereHasTerms = (asphere) => Boolean(asphere) && (
    Math.abs(toNumber(asphere.k ?? asphere.conicK) || 0) > 1e-12
    || Math.abs(toNumber(asphere.A4) || 0) > 1e-18
    || Math.abs(toNumber(asphere.A6) || 0) > 1e-21
    || Math.abs(toNumber(asphere.A8) || 0) > 1e-24
    || Math.abs(toNumber(asphere.A10) || 0) > 1e-27
  );
  const surfaceFeatureFlags = (surfaces = []) => {
    const hasAsphere = surfaces.some((surface) => surface.asphere?.active === true || asphereHasTerms(surface.asphere));
    const hasTilt = surfaces.some((surface) => (
      Math.abs(toNumber(surface.tiltY) || 0) > 1e-9
      || Math.abs(toNumber(surface.tiltZ) || 0) > 1e-9
    ));
    return {
      hasAsphere,
      hasTilt,
      unsupportedByWorker: hasAsphere || hasTilt,
      unsupportedFeatures: [
        hasAsphere ? "active asphere" : "",
        hasTilt ? "surface tilt" : ""
      ].filter(Boolean)
    };
  };
  const surfaceSignature = (surfaces = []) => JSON.stringify(surfaces.map((surface) => ({
    x: Number((toNumber(surface.x) || 0).toFixed(7)),
    radius: Number((isPlanoRadius(surface.radius) ? 0 : toNumber(surface.radius) || 0).toFixed(7)),
    semiDiameter: Number((toNumber(surface.semiDiameter) || 0).toFixed(7)),
    nBefore: Number((toNumber(surface.nBefore) || 1).toFixed(8)),
    nAfter: Number((toNumber(surface.nAfter) || 1).toFixed(8)),
    isStop: surface.isStop === true,
    decenterY: Number((toNumber(surface.decenterY) || 0).toFixed(7)),
    decenterZ: Number((toNumber(surface.decenterZ) || 0).toFixed(7)),
    tiltY: Number((toNumber(surface.tiltY) || 0).toFixed(7)),
    tiltZ: Number((toNumber(surface.tiltZ) || 0).toFixed(7)),
    asphere: surface.asphere ? {
      k: Number((toNumber(surface.asphere.k ?? surface.asphere.conicK) || 0).toPrecision(12)),
      A4: Number((toNumber(surface.asphere.A4) || 0).toPrecision(12)),
      A6: Number((toNumber(surface.asphere.A6) || 0).toPrecision(12)),
      A8: Number((toNumber(surface.asphere.A8) || 0).toPrecision(12)),
      A10: Number((toNumber(surface.asphere.A10) || 0).toPrecision(12))
    } : null
  })));

  const fftRadix2 = (realInput, imagInput = null, inverse = false) => {
    const n = realInput.length;
    const real = realInput.slice();
    const imag = imagInput ? imagInput.slice() : new Array(n).fill(0);
    for (let i = 1, j = 0; i < n; i += 1) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
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
        for (let j = 0; j < len / 2; j += 1) {
          const uReal = real[i + j];
          const uImag = imag[i + j];
          const vReal = real[i + j + len / 2] * wReal - imag[i + j + len / 2] * wImag;
          const vImag = real[i + j + len / 2] * wImag + imag[i + j + len / 2] * wReal;
          real[i + j] = uReal + vReal;
          imag[i + j] = uImag + vImag;
          real[i + j + len / 2] = uReal - vReal;
          imag[i + j + len / 2] = uImag - vImag;
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
    return { real, imag };
  };

  const fieldDirection3D = (fieldAngleDegrees = 0, orientation = "tangential") => {
    const angle = (toNumber(fieldAngleDegrees) || 0) * Math.PI / 180;
    const axial = -Math.cos(angle);
    const transverse = Math.sin(angle);
    if (orientation === "sagittal") return normalizeVector3(axial, 0, transverse) || { x: -1, y: 0, z: 0 };
    if (orientation === "diagonal") {
      const component = transverse / Math.SQRT2;
      return normalizeVector3(axial, component, component) || { x: -1, y: 0, z: 0 };
    }
    return normalizeVector3(axial, transverse, 0) || { x: -1, y: 0, z: 0 };
  };

  const rayTrace3DAxes = (orientation = "tangential") => {
    if (orientation === "sagittal") return { tangential: { y: 0, z: 1 }, sagittal: { y: 1, z: 0 } };
    if (orientation === "diagonal") {
      return {
        tangential: { y: 1 / Math.SQRT2, z: 1 / Math.SQRT2 },
        sagittal: { y: -1 / Math.SQRT2, z: 1 / Math.SQRT2 }
      };
    }
    return { tangential: { y: 1, z: 0 }, sagittal: { y: 0, z: 1 } };
  };
  const project = (point, axis) => point.y * axis.y + point.z * axis.z;

  const apertureMiss3D = (point, surface) => {
    const dy = point.y - (toNumber(surface.decenterY) || 0);
    const dz = point.z - (toNumber(surface.decenterZ) || 0);
    return Math.hypot(dy, dz) > (toNumber(surface.semiDiameter) || 0) + 1e-7;
  };

  const intersectRayWithSurface3D = (ray, surface) => {
    const radius = toNumber(surface.radius);
    if (isPlanoRadius(radius)) {
      if (Math.abs(ray.dx) < 1e-9) return { status: "invalid" };
      const t = (surface.x - ray.x) / ray.dx;
      if (!(t > 1e-7)) return { status: "invalid" };
      const hitPoint = { x: surface.x, y: ray.y + ray.dy * t, z: ray.z + ray.dz * t };
      if (apertureMiss3D(hitPoint, surface)) return { status: "missed aperture", hitPoint };
      return { status: "valid", hitPoint, t };
    }
    const centerX = surface.x - radius;
    const centerY = toNumber(surface.decenterY) || 0;
    const centerZ = toNumber(surface.decenterZ) || 0;
    const ox = ray.x - centerX;
    const oy = ray.y - centerY;
    const oz = ray.z - centerZ;
    const b = 2 * (ox * ray.dx + oy * ray.dy + oz * ray.dz);
    const c = ox ** 2 + oy ** 2 + oz ** 2 - radius ** 2;
    const disc = b ** 2 - 4 * c;
    if (disc < 0) return { status: "invalid" };
    const sqrtDisc = Math.sqrt(Math.max(0, disc));
    const candidates = [(-b - sqrtDisc) / 2, (-b + sqrtDisc) / 2]
      .filter((candidate) => candidate > 1e-7)
      .sort((a, bValue) => a - bValue);
    if (!candidates.length) return { status: "invalid" };
    const t = candidates[0];
    const hitPoint = { x: ray.x + ray.dx * t, y: ray.y + ray.dy * t, z: ray.z + ray.dz * t };
    if (apertureMiss3D(hitPoint, surface)) return { status: "missed aperture", hitPoint };
    return { status: "valid", hitPoint, t };
  };

  const surfaceNormal3D = (ray, surface, hitPoint) => {
    const radius = toNumber(surface.radius);
    let normal;
    if (isPlanoRadius(radius)) {
      normal = { x: ray.dx < 0 ? -1 : 1, y: 0, z: 0 };
    } else {
      normal = normalizeVector3(
        hitPoint.x - (surface.x - radius),
        hitPoint.y - (toNumber(surface.decenterY) || 0),
        hitPoint.z - (toNumber(surface.decenterZ) || 0)
      );
    }
    if (!normal) return null;
    if (dot3({ x: ray.dx, y: ray.dy, z: ray.dz }, normal) < 0) {
      normal = { x: -normal.x, y: -normal.y, z: -normal.z };
    }
    return normal;
  };

  const refractRay3D = (ray, surface, hitPoint) => {
    const normal = surfaceNormal3D(ray, surface, hitPoint);
    if (!normal) return { status: "invalid" };
    const incoming = { x: ray.dx, y: ray.dy, z: ray.dz };
    const cosIncoming = dot3(incoming, normal);
    const tangent = {
      x: incoming.x - cosIncoming * normal.x,
      y: incoming.y - cosIncoming * normal.y,
      z: incoming.z - cosIncoming * normal.z
    };
    const nBefore = toNumber(surface.nBefore) || 1;
    const nAfter = toNumber(surface.nAfter) || 1;
    const eta = nBefore / nAfter;
    const transmittedTangent = { x: tangent.x * eta, y: tangent.y * eta, z: tangent.z * eta };
    const normalTermSquared = 1 - (transmittedTangent.x ** 2 + transmittedTangent.y ** 2 + transmittedTangent.z ** 2);
    if (normalTermSquared < -1e-7) return { status: "total internal reflection" };
    const normalTerm = Math.sqrt(Math.max(0, normalTermSquared));
    const direction = normalizeVector3(
      transmittedTangent.x + normal.x * normalTerm,
      transmittedTangent.y + normal.y * normalTerm,
      transmittedTangent.z + normal.z * normalTerm
    );
    if (!direction) return { status: "invalid" };
    return { status: "valid", ray: { x: hitPoint.x, y: hitPoint.y, z: hitPoint.z, dx: direction.x, dy: direction.y, dz: direction.z } };
  };

  const traceRay3D = (ray, surfaces) => {
    const direction = normalizeVector3(ray?.dx, ray?.dy, ray?.dz);
    if (!direction) return { status: "invalid", finalRay: ray, path: [] };
    let currentRay = { ...ray, dx: direction.x, dy: direction.y, dz: direction.z };
    const path = [{ x: currentRay.x, y: currentRay.y, z: currentRay.z }];
    for (const surface of surfaces) {
      const intersection = intersectRayWithSurface3D(currentRay, surface);
      if (intersection.status !== "valid") {
        return { status: intersection.status, finalRay: currentRay, failedSurface: surface, hitPoint: intersection.hitPoint, path };
      }
      path.push(intersection.hitPoint);
      if (surface.isStop || Math.abs((toNumber(surface.nBefore) || 1) - (toNumber(surface.nAfter) || 1)) < 1e-7) {
        currentRay = { ...intersection.hitPoint, dx: currentRay.dx, dy: currentRay.dy, dz: currentRay.dz };
        continue;
      }
      const refracted = refractRay3D(currentRay, surface, intersection.hitPoint);
      if (refracted.status !== "valid") {
        return { status: refracted.status, finalRay: currentRay, failedSurface: surface, hitPoint: intersection.hitPoint, path };
      }
      currentRay = refracted.ray;
    }
    return { status: "valid", finalRay: currentRay, path };
  };

  const imagePlaneIntersection3D = (ray, imagePlaneX = 0) => {
    if (!ray || Math.abs(ray.dx) < 1e-9) return null;
    const t = (imagePlaneX - ray.x) / ray.dx;
    if (!(t >= 0)) return null;
    return { x: imagePlaneX, y: ray.y + ray.dy * t, z: ray.z + ray.dz * t };
  };

  const makeRayFromPupil = ({ sample, fieldAngleDegrees, orientation, apertureRadius, referenceX, startX }) => {
    const direction = fieldDirection3D(fieldAngleDegrees, orientation);
    const apertureY = sample.pupilU * apertureRadius;
    const apertureZ = sample.pupilV * apertureRadius;
    const tToReference = Math.abs(direction.x) > 1e-9 ? (referenceX - startX) / direction.x : 0;
    return {
      x: startX,
      y: apertureY - direction.y * tToReference,
      z: apertureZ - direction.z * tToReference,
      dx: direction.x,
      dy: direction.y,
      dz: direction.z,
      apertureY,
      apertureZ,
      pupilU: sample.pupilU,
      pupilV: sample.pupilV,
      pupilWeight: sample.pupilWeight,
      sourceWeight: sample.sourceWeight,
      isChiefReference: sample.isChiefReference === true
    };
  };

  const traceChiefImageHeight = (surfaces, setup, fieldAngleDegrees) => {
    const ray = makeRayFromPupil({
      sample: { pupilU: 0, pupilV: 0, pupilWeight: 0, sourceWeight: 0, isChiefReference: true },
      fieldAngleDegrees,
      orientation: setup.orientation,
      apertureRadius: setup.apertureRadius,
      referenceX: setup.referenceX,
      startX: setup.startX
    });
    const traced = traceRay3D(ray, surfaces);
    const imagePoint = traced.status === "valid" ? imagePlaneIntersection3D(traced.finalRay, setup.imagePlaneX) : null;
    return {
      status: imagePoint ? "valid" : traced.status,
      imagePoint,
      imageHeight: imagePoint ? Math.hypot(imagePoint.y, imagePoint.z) : NaN
    };
  };

  const solveFieldAngleForImageHeight = (surfaces, setup, targetImageHeight) => {
    if (targetImageHeight <= 1e-6) return { status: "solved", fieldAngleDegrees: 0, imageHeight: 0, residualMm: 0 };
    const maxAngle = Math.max(1, Math.min(70, toNumber(setup.maxFieldAngleDegrees) || 45));
    const scanned = [];
    for (let angle = 0; angle <= maxAngle + 1e-9; angle += 1) {
      const sample = traceChiefImageHeight(surfaces, setup, angle);
      if (sample.status === "valid" && Number.isFinite(sample.imageHeight)) scanned.push({ angle, ...sample });
    }
    if (scanned.length < 2) return { status: "invalid", fieldAngleDegrees: NaN, imageHeight: NaN, residualMm: NaN };
    let bracket = null;
    for (let i = 1; i < scanned.length; i += 1) {
      const low = scanned[i - 1];
      const high = scanned[i];
      if (high.imageHeight + 1e-6 < low.imageHeight) continue;
      if (low.imageHeight <= targetImageHeight && targetImageHeight <= high.imageHeight) {
        bracket = { low, high };
        break;
      }
    }
    if (!bracket) {
      const best = scanned.reduce((selected, sample) => (
        Math.abs(sample.imageHeight - targetImageHeight) < Math.abs(selected.imageHeight - targetImageHeight) ? sample : selected
      ), scanned[0]);
      return { status: "unreachable", fieldAngleDegrees: best.angle, imageHeight: best.imageHeight, residualMm: Math.abs(best.imageHeight - targetImageHeight) };
    }
    let low = bracket.low.angle;
    let high = bracket.high.angle;
    for (let i = 0; i < 18; i += 1) {
      const mid = (low + high) / 2;
      const sample = traceChiefImageHeight(surfaces, setup, mid);
      if (sample.status === "valid" && sample.imageHeight >= targetImageHeight) high = mid;
      else low = mid;
    }
    const final = traceChiefImageHeight(surfaces, setup, high);
    return {
      status: final.status === "valid" ? "solved" : "invalid",
      fieldAngleDegrees: high,
      imageHeight: final.imageHeight,
      residualMm: Number.isFinite(final.imageHeight) ? Math.abs(final.imageHeight - targetImageHeight) : NaN
    };
  };

  const createBaseCells = (gridSize) => {
    const cellWidth = 2 / gridSize;
    const cells = [];
    for (let row = 0; row < gridSize; row += 1) {
      for (let col = 0; col < gridSize; col += 1) {
        const pupilU = -1 + (col + 0.5) * cellWidth;
        const pupilV = -1 + (row + 0.5) * cellWidth;
        if (Math.hypot(pupilU, pupilV) <= 1) {
          cells.push({ row, col, pupilU, pupilV, width: cellWidth, sourceWeight: cellWidth * cellWidth / Math.PI });
        }
      }
    }
    return cells;
  };

  const classifyCell = (surfaces, setup, cell, fieldAngleDegrees) => {
    const ray = makeRayFromPupil({
      sample: cell,
      fieldAngleDegrees,
      orientation: setup.orientation,
      apertureRadius: setup.apertureRadius,
      referenceX: setup.referenceX,
      startX: setup.startX
    });
    const traced = traceRay3D(ray, surfaces);
    const imagePoint = traced.status === "valid" ? imagePlaneIntersection3D(traced.finalRay, setup.imagePlaneX) : null;
    const status = traced.status === "valid" && !imagePoint ? "invalid" : traced.status;
    return { ...cell, ray, traced, imagePoint, status };
  };

  const adaptivePupilSamples = (surfaces, setup, quality, fieldAngleDegrees) => {
    const profile = QUALITY_PROFILES[quality] || QUALITY_PROFILES.interactive;
    const baseCells = createBaseCells(profile.baseGrid);
    const classified = baseCells.map((cell) => classifyCell(surfaces, setup, cell, fieldAngleDegrees));
    const byKey = new Map(classified.map((cell) => [`${cell.row}:${cell.col}`, cell.status]));
    const boundaryCells = classified.filter((cell) => {
      const neighbors = [[1, 0], [-1, 0], [0, 1], [0, -1]]
        .map(([dr, dc]) => byKey.get(`${cell.row + dr}:${cell.col + dc}`))
        .filter(Boolean);
      return neighbors.some((status) => status !== cell.status) || cell.status !== "valid";
    });
    const boundarySet = new Set(boundaryCells.map((cell) => `${cell.row}:${cell.col}`));
    const refined = [];
    classified.forEach((cell) => {
      if (!boundarySet.has(`${cell.row}:${cell.col}`)) {
        refined.push(cell);
        return;
      }
      const childWidth = cell.width / 2;
      [-0.25, 0.25].forEach((du) => {
        [-0.25, 0.25].forEach((dv) => {
          const child = {
            row: cell.row,
            col: cell.col,
            pupilU: cell.pupilU + du * cell.width,
            pupilV: cell.pupilV + dv * cell.width,
            width: childWidth,
            sourceWeight: cell.sourceWeight / 4,
            refinedFromBoundary: true
          };
          if (Math.hypot(child.pupilU, child.pupilV) <= 1) refined.push(classifyCell(surfaces, setup, child, fieldAngleDegrees));
        });
      });
    });
    const sourcePupilWeight = refined.reduce((sum, sample) => sum + sample.sourceWeight, 0) || 1;
    const valid = refined.filter((sample) => sample.status === "valid" && sample.imagePoint);
    const survivingPupilWeight = valid.reduce((sum, sample) => sum + sample.sourceWeight, 0);
    const clippedPupilWeight = refined
      .filter((sample) => sample.status !== "valid")
      .reduce((sum, sample) => sum + sample.sourceWeight, 0);
    const normalizedDenominator = survivingPupilWeight || 1;
    const validWithWeights = valid.map((sample) => ({ ...sample, pupilWeight: sample.sourceWeight / normalizedDenominator }));
    const effectiveSampleCount = validWithWeights.length
      ? 1 / validWithWeights.reduce((sum, sample) => sum + sample.pupilWeight ** 2, 0)
      : 0;
    return {
      validSamples: validWithWeights,
      allSamples: refined,
      diagnostics: {
        sourcePupilWeight,
        survivingPupilWeight,
        clippedPupilWeight,
        transmittedFraction: sourcePupilWeight > 0 ? survivingPupilWeight / sourcePupilWeight : 0,
        effectiveSampleCount,
        boundaryRefinementCount: refined.filter((sample) => sample.refinedFromBoundary).length,
        validEnergyRayCount: validWithWeights.length,
        clippedEnergyRayCount: refined.length - validWithWeights.length
      }
    };
  };

  const sharedLsfGrid = (qualitySamples, maxFrequencyLpMm) => {
    const allCoordinates = qualitySamples.flatMap((item) => item.samples.map((sample) => sample.coordinate)).filter(Number.isFinite);
    if (allCoordinates.length < 2) return null;
    const min = Math.min(...allCoordinates);
    const max = Math.max(...allCoordinates);
    const span = Math.max(max - min, 1 / Math.max(1, maxFrequencyLpMm));
    const padding = Math.max(span * 0.2, 1 / Math.max(1, maxFrequencyLpMm));
    const left = min - padding;
    const right = max + padding;
    const binPitchMm = 1 / (Math.max(1, maxFrequencyLpMm) * 16);
    const binCount = clamp(nextPowerOfTwo(Math.ceil((right - left) / binPitchMm) + 1), 256, 16384);
    const actualPitch = (right - left) / (binCount - 1);
    const fftSize = nextPowerOfTwo(binCount * 4);
    return {
      left,
      right,
      binPitchMm: actualPitch,
      binCount,
      fftSize,
      nyquistFrequencyLpMm: 1 / (2 * actualPitch),
      samplesPerCycleAtMax: 1 / (Math.max(1, maxFrequencyLpMm) * actualPitch),
      sharedLsfGridUsed: true
    };
  };

  const mtfFromSamplesOnGrid = (samples, grid, frequencies) => {
    const bins = new Array(grid.binCount).fill(0);
    samples.forEach((sample) => {
      if (!Number.isFinite(sample.coordinate) || !(sample.weight > 0)) return;
      const position = (sample.coordinate - grid.left) / grid.binPitchMm;
      const index = Math.floor(position);
      const fraction = position - index;
      if (index >= 0 && index < bins.length) bins[index] += sample.weight * (1 - fraction);
      if (index + 1 >= 0 && index + 1 < bins.length) bins[index + 1] += sample.weight * fraction;
    });
    const energy = bins.reduce((sum, value) => sum + value, 0) || 1;
    const normalized = bins.map((value) => value / energy);
    const padded = normalized.concat(new Array(grid.fftSize - normalized.length).fill(0));
    const fft = fftRadix2(padded);
    const dc = Math.hypot(fft.real[0], fft.imag[0]) || 1;
    const values = {};
    const diagnostics = {};
    frequencies.forEach((frequency) => {
      const bin = frequency * grid.fftSize * grid.binPitchMm;
      const leftIndex = Math.floor(bin);
      const fraction = bin - leftIndex;
      let value = NaN;
      if (leftIndex >= 0 && leftIndex + 1 < fft.real.length) {
        const a = Math.hypot(fft.real[leftIndex], fft.imag[leftIndex]) / dc;
        const b = Math.hypot(fft.real[leftIndex + 1], fft.imag[leftIndex + 1]) / dc;
        value = a + (b - a) * fraction;
        const kernelMtf = sinc(frequency * grid.binPitchMm) ** 2;
        if (frequency < 0.5 * grid.nyquistFrequencyLpMm && kernelMtf > 0.25) {
          value = clamp(value / kernelMtf, 0, 1);
        }
        diagnostics[frequency] = {
          kernelMtf,
          correctionApplied: frequency < 0.5 * grid.nyquistFrequencyLpMm && kernelMtf > 0.25,
          frequencyLimited: frequency > 0.8 * grid.nyquistFrequencyLpMm
        };
      }
      values[frequency] = value;
    });
    return { values, diagnostics };
  };

  const mtfValueAtFrequency = (curve, frequency) => {
    const target = toNumber(frequency);
    const frequencies = curve?.frequencies || [];
    const values = curve?.values || [];
    if (!Number.isFinite(target) || !frequencies.length || frequencies.length !== values.length) return NaN;
    if (target <= frequencies[0]) return values[0];
    for (let index = 1; index < frequencies.length; index += 1) {
      if (target <= frequencies[index]) {
        const leftFrequency = frequencies[index - 1];
        const rightFrequency = frequencies[index];
        const leftValue = values[index - 1];
        const rightValue = values[index];
        const span = rightFrequency - leftFrequency || 1;
        return leftValue + (rightValue - leftValue) * ((target - leftFrequency) / span);
      }
    }
    return values[values.length - 1];
  };

  const mainPanelQuality = (quality) => MAIN_PANEL_QUALITY_PROFILES[quality] ? quality : "interactive";

  const weightedEqualAreaPupilSamples = (quality = "interactive", apertureRadius = 10) => {
    const profile = MAIN_PANEL_QUALITY_PROFILES[mainPanelQuality(quality)];
    const rings = profile.rings;
    const radius = Math.max(0.1, toNumber(apertureRadius) || 10);
    const samples = [{
      y: 0,
      z: 0,
      pupilU: 0,
      pupilV: 0,
      ring: 0,
      isChiefReference: true,
      weight: 0,
      pupilWeight: 0
    }];
    for (let ring = 1; ring <= rings; ring += 1) {
      const innerRadius = Math.sqrt((ring - 1) / rings);
      const outerRadius = Math.sqrt(ring / rings);
      const normalizedRadius = Math.sqrt((innerRadius ** 2 + outerRadius ** 2) / 2);
      const angularCount = Math.max(8, ring * 6);
      const annulusArea = Math.PI * (outerRadius ** 2 - innerRadius ** 2);
      const sampleWeight = annulusArea / angularCount;
      for (let index = 0; index < angularCount; index += 1) {
        const angle = 2 * Math.PI * (index + 0.5) / angularCount;
        const pupilU = normalizedRadius * Math.cos(angle);
        const pupilV = normalizedRadius * Math.sin(angle);
        samples.push({
          y: radius * pupilU,
          z: radius * pupilV,
          pupilU,
          pupilV,
          ring,
          innerRadius,
          outerRadius,
          annulusArea,
          isChiefReference: false,
          weight: sampleWeight,
          pupilWeight: sampleWeight
        });
      }
    }
    const energyWeight = samples
      .filter((sample) => !sample.isChiefReference)
      .reduce((sum, sample) => sum + sample.pupilWeight, 0) || 1;
    return samples.map((sample) => {
      const pupilWeight = sample.isChiefReference ? 0 : sample.pupilWeight / energyWeight;
      return { ...sample, weight: pupilWeight, pupilWeight };
    });
  };

  const makeWeightedRayBundle = (setup, fieldAngleDegrees = 0, quality = "interactive") => {
    const direction = fieldDirection3D(fieldAngleDegrees, setup.orientation || "tangential");
    const apertureRadius = Math.max(0.1, toNumber(setup.apertureRadius) || 10);
    const referenceX = toNumber(setup.referenceX) || 0;
    const startX = toNumber(setup.startX) || referenceX + 25;
    const tToReference = Math.abs(direction.x) > 1e-9 ? (referenceX - startX) / direction.x : 0;
    return weightedEqualAreaPupilSamples(quality, apertureRadius).map((sample) => ({
      x: startX,
      y: sample.y - direction.y * tToReference,
      z: sample.z - direction.z * tToReference,
      dx: direction.x,
      dy: direction.y,
      dz: direction.z,
      apertureY: sample.y,
      apertureZ: sample.z,
      pupilU: sample.pupilU,
      pupilV: sample.pupilV,
      pupilWeight: sample.pupilWeight,
      pupilRing: sample.ring,
      isChiefReference: sample.isChiefReference === true
    }));
  };

  const traceWeightedPupilForMainMtf = (surfaces, setup, fieldAngleDegrees, quality = "interactive") => {
    const rays = makeWeightedRayBundle(setup, fieldAngleDegrees, quality).map((ray) => {
      const traced = traceRay3D(ray, surfaces);
      const imagePoint = traced.status === "valid" ? imagePlaneIntersection3D(traced.finalRay, setup.imagePlaneX) : null;
      const status = traced.status === "valid" && !imagePoint ? "invalid" : traced.status;
      return {
        ...traced,
        status,
        inputRay: ray,
        imagePoint,
        weight: ray.pupilWeight,
        isChiefReference: ray.isChiefReference === true
      };
    });
    const validRays = rays.filter((ray) => ray.status === "valid" && ray.imagePoint);
    const energyRays = rays.filter((ray) => ray.weight > 0);
    const validEnergyRays = validRays.filter((ray) => ray.weight > 0);
    return {
      rays,
      validRays,
      totalRayCount: rays.length,
      totalEnergyRayCount: energyRays.length,
      validEnergyRayCount: validEnergyRays.length,
      chiefReferenceRayCount: rays.filter((ray) => ray.isChiefReference).length,
      clippedRayCount: rays.filter((ray) => ray.status === "missed aperture").length,
      clippedEnergyRayCount: energyRays.filter((ray) => ray.status === "missed aperture").length
    };
  };

  const buildWeightedLsf = (samples, maxFrequencyLpMm) => {
    const energySamples = samples.filter((sample) => (
      Number.isFinite(sample.coordinate)
      && Number.isFinite(sample.weight)
      && sample.weight > 0
    ));
    const coordinates = energySamples.map((sample) => sample.coordinate);
    if (coordinates.length < 2) return { status: "invalid", warning: "Too few valid ray intercepts for LSF.", bins: [], binPitchMm: NaN };
    const maxFrequency = Math.max(10, toNumber(maxFrequencyLpMm) || 100);
    const samplesPerCycleTarget = 8;
    const minPitch = 1 / (Math.max(2, maxFrequency) * samplesPerCycleTarget);
    const minCoordinate = Math.min(...coordinates);
    const maxCoordinate = Math.max(...coordinates);
    const span = Math.max(maxCoordinate - minCoordinate, minPitch * 16);
    const padding = Math.max(span * 0.25, minPitch * 8);
    const left = minCoordinate - padding;
    const right = maxCoordinate + padding;
    const rawBinCount = Math.ceil((right - left) / minPitch) + 1;
    const maxBinCount = 8192;
    const requestedBinCount = nextPowerOfTwo(rawBinCount);
    const binCount = clamp(requestedBinCount, 128, maxBinCount);
    const binPitchMm = (right - left) / (binCount - 1);
    const bins = new Array(binCount).fill(0);
    energySamples.forEach((sample) => {
      const position = (sample.coordinate - left) / binPitchMm;
      const index = Math.floor(position);
      const fraction = position - index;
      if (index >= 0 && index < binCount) bins[index] += sample.weight * (1 - fraction);
      if (index + 1 >= 0 && index + 1 < binCount) bins[index + 1] += sample.weight * fraction;
    });
    const energy = bins.reduce((sum, value) => sum + value, 0);
    if (!(energy > 0)) return { status: "invalid", warning: "LSF energy is zero.", bins, binPitchMm };
    const nyquistFrequencyLpMm = 1 / (2 * binPitchMm);
    const samplesPerCycleAtMax = 1 / (maxFrequency * binPitchMm);
    const reliableMaxFrequencyLpMm = 0.8 * nyquistFrequencyLpMm;
    const frequencyLimited = requestedBinCount > maxBinCount || maxFrequency > reliableMaxFrequencyLpMm;
    return {
      status: "valid",
      bins: bins.map((value) => value / energy),
      binPitchMm,
      supportMm: right - left,
      binCount,
      requestedBinCount,
      maxBinCount,
      nyquistFrequencyLpMm,
      samplesPerCycleAtMax,
      samplesPerCycleTarget,
      reliableMaxFrequencyLpMm,
      effectiveMaxFrequencyLpMm: frequencyLimited ? reliableMaxFrequencyLpMm : maxFrequency,
      frequencyLimited,
      energySampleCount: energySamples.length
    };
  };

  const axisMtfResultFromWeightedLsf = (axis, samples, options = {}) => {
    const maxFrequency = Math.max(10, toNumber(options.maxFrequencyLpMm) || 100);
    const frequencyStep = Math.max(1, toNumber(options.frequencyStepLpMm) || 5);
    const lsf = buildWeightedLsf(samples, maxFrequency);
    if (lsf.status !== "valid") {
      return {
        axis,
        engine: "geometricLsfFft",
        status: lsf.status,
        warning: lsf.warning,
        frequencies: [],
        values: [],
        readouts: Object.fromEntries(MTF_READOUT_FREQUENCIES.map((frequency) => [frequency, NaN])),
        lsf
      };
    }
    const paddedLength = nextPowerOfTwo(lsf.bins.length * 2);
    const padded = lsf.bins.concat(new Array(paddedLength - lsf.bins.length).fill(0));
    const fft = fftRadix2(padded);
    const dc = Math.hypot(fft.real[0], fft.imag[0]) || 1;
    const effectiveMaxFrequency = Math.max(0, Math.min(maxFrequency, lsf.effectiveMaxFrequencyLpMm || maxFrequency));
    const frequencies = [];
    const values = [];
    const frequencyDiagnostics = {};
    for (let frequency = 0; frequency <= effectiveMaxFrequency + 1e-6; frequency += frequencyStep) {
      const roundedFrequency = Number(frequency.toFixed(4));
      const bin = roundedFrequency * paddedLength * lsf.binPitchMm;
      const leftIndex = Math.floor(bin);
      const rightIndex = Math.min(fft.real.length - 1, leftIndex + 1);
      const t = bin - leftIndex;
      const magnitudeAt = (index) => Math.hypot(fft.real[index] || 0, fft.imag[index] || 0) / dc;
      const rawValue = leftIndex >= 0 && leftIndex < fft.real.length
        ? magnitudeAt(leftIndex) * (1 - t) + magnitudeAt(rightIndex) * t
        : NaN;
      const kernelMtf = sinc(roundedFrequency * lsf.binPitchMm) ** 2;
      const correctionApplied = roundedFrequency < 0.5 * lsf.nyquistFrequencyLpMm && kernelMtf > 0.25;
      const value = correctionApplied ? rawValue / kernelMtf : rawValue;
      frequencies.push(roundedFrequency);
      values.push(clamp(value, 0, 1));
      frequencyDiagnostics[roundedFrequency] = {
        rawValue,
        kernelMtf,
        correctionApplied,
        frequencyLimited: roundedFrequency > lsf.reliableMaxFrequencyLpMm
      };
    }
    const result = {
      axis,
      engine: "geometricLsfFft",
      status: "valid",
      frequencies,
      values,
      lsf: {
        binPitchMm: lsf.binPitchMm,
        binCount: lsf.binCount,
        supportMm: lsf.supportMm,
        kernel: "linear B-spline / triangular binning",
        kernelCorrection: "sinc(f * binPitch)^2 correction below 50% Nyquist when stable",
        frequencyDiagnostics,
        nyquistFrequencyLpMm: lsf.nyquistFrequencyLpMm,
        samplesPerCycleAtMax: lsf.samplesPerCycleAtMax,
        samplesPerCycleTarget: lsf.samplesPerCycleTarget,
        zeroPaddingFactor: paddedLength / lsf.binCount,
        energySampleCount: lsf.energySampleCount,
        requestedMaxFrequencyLpMm: maxFrequency,
        effectiveMaxFrequencyLpMm: effectiveMaxFrequency,
        frequencyLimited: lsf.frequencyLimited === true,
        warning: lsf.frequencyLimited ? "Frequency range limited by LSF sampling." : ""
      },
      warnings: lsf.frequencyLimited ? ["Frequency range limited by LSF sampling."] : []
    };
    return {
      ...result,
      readouts: Object.fromEntries(MTF_READOUT_FREQUENCIES.map((frequency) => [
        frequency,
        mtfValueAtFrequency(result, frequency)
      ]))
    };
  };

  const calculateMainGeometricLsfMtf = (surfaces, setup, request = {}) => {
    const quality = mainPanelQuality(request.quality || setup.quality);
    const fieldAngleDegrees = toNumber(request.fieldAngleDegrees) || 0;
    const traced = traceWeightedPupilForMainMtf(surfaces, setup, fieldAngleDegrees, quality);
    const warnings = [];
    if (traced.validEnergyRayCount < 12) warnings.push(`Fewer than 12 valid weighted pupil rays for ${fieldAngleDegrees}° d-line LSF/FFT MTF.`);
    if (traced.totalEnergyRayCount && traced.clippedEnergyRayCount / traced.totalEnergyRayCount > 0.5) warnings.push(`More than 50% of weighted pupil rays are clipped at ${fieldAngleDegrees}° d-line.`);
    const chiefRay = traced.validRays.reduce((closest, ray) => {
      const pupilDistance = Math.hypot(ray.inputRay?.pupilU || 0, ray.inputRay?.pupilV || 0);
      const closestDistance = closest ? Math.hypot(closest.inputRay?.pupilU || 0, closest.inputRay?.pupilV || 0) : Infinity;
      return pupilDistance < closestDistance ? ray : closest;
    }, null);
    const chiefY = chiefRay?.imagePoint?.y ?? NaN;
    const chiefZ = chiefRay?.imagePoint?.z ?? NaN;
    if (traced.validEnergyRayCount < 4 || !chiefRay) {
      return {
        fieldKey: request.fieldKey,
        fieldName: request.fieldName || "Field",
        fieldAngleDegrees,
        spectralLineKey: request.spectralLineKey || "d",
        wavelengthNm: toNumber(request.wavelengthNm) || 587.6,
        imagePlaneX: setup.imagePlaneX,
        planeLabel: request.planeLabel || "Current sensor plane",
        status: "too-few-rays",
        engine: "geometricLsfFft",
        engineLabel: "Geometric LSF/FFT preview — not physical wavefront MTF",
        quality,
        isFallback: false,
        fallbackReason: "",
        validRayCount: traced.validRays.length,
        totalRayCount: traced.totalRayCount,
        validEnergyRayCount: traced.validEnergyRayCount,
        totalEnergyRayCount: traced.totalEnergyRayCount,
        chiefReferenceRayCount: traced.chiefReferenceRayCount,
        clippedRayCount: traced.clippedRayCount,
        clippedEnergyRayCount: traced.clippedEnergyRayCount,
        apertureStopX: setup.referenceX,
        chiefY,
        chiefZ,
        tangential: axisMtfResultFromWeightedLsf("tangential", [], request),
        sagittal: axisMtfResultFromWeightedLsf("sagittal", [], request),
        combinedRms: NaN,
        lsfKernel: "linear B-spline / triangular binning",
        warnings
      };
    }
    const axes = rayTrace3DAxes(setup.orientation || "tangential");
    const weightedSamples = traced.validRays.map((ray) => {
      const dy = ray.imagePoint.y - chiefY;
      const dz = ray.imagePoint.z - chiefZ;
      return {
        tangential: project({ y: dy, z: dz }, axes.tangential),
        sagittal: project({ y: dy, z: dz }, axes.sagittal),
        weight: Number.isFinite(ray.weight) ? ray.weight : 1
      };
    });
    const rmsFor = (axis) => {
      const totalWeight = weightedSamples.reduce((sum, sample) => sum + sample.weight, 0) || 1;
      return Math.sqrt(weightedSamples.reduce((sum, sample) => sum + sample.weight * sample[axis] ** 2, 0) / totalWeight);
    };
    const tangential = axisMtfResultFromWeightedLsf("tangential", weightedSamples.map((sample) => ({
      coordinate: sample.tangential,
      weight: sample.weight
    })), request);
    const sagittal = axisMtfResultFromWeightedLsf("sagittal", weightedSamples.map((sample) => ({
      coordinate: sample.sagittal,
      weight: sample.weight
    })), request);
    const tangentialSigma = rmsFor("tangential");
    const sagittalSigma = rmsFor("sagittal");
    const combinedRms = Number.isFinite(tangentialSigma) && Number.isFinite(sagittalSigma)
      ? Math.sqrt((tangentialSigma ** 2 + sagittalSigma ** 2) / 2)
      : NaN;
    return {
      fieldKey: request.fieldKey,
      fieldName: request.fieldName || "Field",
      fieldAngleDegrees,
      spectralLineKey: request.spectralLineKey || "d",
      wavelengthNm: toNumber(request.wavelengthNm) || 587.6,
      imagePlaneX: setup.imagePlaneX,
      planeLabel: request.planeLabel || "Current sensor plane",
      status: traced.validRays.length >= 12 && tangential.status === "valid" && sagittal.status === "valid" ? "valid" : "too-few-rays",
      engine: "geometricLsfFft",
      engineLabel: "Geometric LSF/FFT preview — not physical wavefront MTF",
      quality,
      isFallback: false,
      fallbackReason: "",
      validRayCount: traced.validRays.length,
      totalRayCount: traced.totalRayCount,
      validEnergyRayCount: traced.validEnergyRayCount,
      totalEnergyRayCount: traced.totalEnergyRayCount,
      chiefReferenceRayCount: traced.chiefReferenceRayCount,
      clippedRayCount: traced.clippedRayCount,
      clippedEnergyRayCount: traced.clippedEnergyRayCount,
      apertureStopX: setup.referenceX,
      chiefY,
      chiefZ,
      tangential: { ...tangential, sigma: tangentialSigma },
      sagittal: { ...sagittal, sigma: sagittalSigma },
      combinedRms,
      lsfKernel: "linear B-spline / triangular binning",
      warnings: [...warnings, ...(tangential.warnings || []), ...(sagittal.warnings || [])].filter((warning, index, list) => warning && list.indexOf(warning) === index)
    };
  };

  const traceQuality = (surfaces, setup, quality, fieldAngleDegrees, chiefImagePoint) => {
    const adaptive = adaptivePupilSamples(surfaces, setup, quality, fieldAngleDegrees);
    const axes = rayTrace3DAxes(setup.orientation);
    const samples = adaptive.validSamples.map((sample) => {
      const dy = sample.imagePoint.y - chiefImagePoint.y;
      const dz = sample.imagePoint.z - chiefImagePoint.z;
      return {
        tangential: project({ y: dy, z: dz }, axes.tangential),
        sagittal: project({ y: dy, z: dz }, axes.sagittal),
        weight: sample.pupilWeight
      };
    });
    return {
      quality,
      diagnostics: adaptive.diagnostics,
      axes: {
        tangential: samples.map((sample) => ({ coordinate: sample.tangential, weight: sample.weight })),
        sagittal: samples.map((sample) => ({ coordinate: sample.sagittal, weight: sample.weight }))
      }
    };
  };

  const convergenceForField = (surfaces, setup, field, frequencies, maxFrequencyLpMm) => {
    const solved = solveFieldAngleForImageHeight(surfaces, setup, field.imageHeightMm);
    if (field.imageHeightMm > 1e-6 && solved.status !== "solved") {
      return {
        field,
        comparisons: [{ fieldKey: field.key, fieldName: field.name, status: solved.status, warning: `Field target ${field.imageHeightMm} mm could not be solved.` }],
        qualityDiagnostics: {},
        status: solved.status
      };
    }
    const fieldAngleDegrees = field.imageHeightMm <= 1e-6 ? 0 : solved.fieldAngleDegrees;
    const chief = traceChiefImageHeight(surfaces, setup, fieldAngleDegrees);
    if (chief.status !== "valid" || !chief.imagePoint) {
      return { field, comparisons: [{ fieldKey: field.key, fieldName: field.name, status: "invalid" }], qualityDiagnostics: {}, status: "invalid" };
    }
    const qualities = ["interactive", "high", "reference"].map((quality) => (
      traceQuality(surfaces, setup, quality, fieldAngleDegrees, chief.imagePoint)
    ));
    const qualityDiagnostics = Object.fromEntries(qualities.map((item) => [item.quality, item.diagnostics]));
    const comparisons = [];
    ["sagittal", "tangential"].forEach((axis) => {
      const axisQualitySamples = qualities.map((item) => ({ quality: item.quality, samples: item.axes[axis] }));
      const grid = sharedLsfGrid(axisQualitySamples, maxFrequencyLpMm);
      if (!grid) {
        comparisons.push({ fieldKey: field.key, fieldName: field.name, axis, status: "invalid", warning: "Shared LSF grid could not be built." });
        return;
      }
      const mtfs = Object.fromEntries(axisQualitySamples.map((item) => [
        item.quality,
        mtfFromSamplesOnGrid(item.samples, grid, frequencies)
      ]));
      frequencies.forEach((frequency) => {
        const interactive = mtfs.interactive.values[frequency];
        const high = mtfs.high.values[frequency];
        const reference = mtfs.reference.values[frequency];
        comparisons.push({
          fieldKey: field.key,
          fieldName: field.name,
          imageHeightMm: field.imageHeightMm,
          fieldAngleDegrees,
          axis,
          frequency,
          interactive,
          high,
          reference,
          interactiveHighDelta: Math.abs(interactive - high),
          highReferenceDelta: Math.abs(high - reference),
          status: [interactive, high, reference].every(Number.isFinite) ? "compared" : "invalid",
          sharedLsfGridUsed: true,
          lsf: {
            binPitchMm: grid.binPitchMm,
            binCount: grid.binCount,
            fftSize: grid.fftSize,
            nyquistFrequencyLpMm: grid.nyquistFrequencyLpMm,
            samplesPerCycleAtMax: grid.samplesPerCycleAtMax,
            kernel: mtfs.reference.diagnostics[frequency]
          }
        });
      });
    });
    return { field, comparisons, qualityDiagnostics, status: "compared" };
  };

  const calculateGeometricLsfConvergence = (payload = {}) => {
    const started = Date.now();
    if (payload.solverContractVersion && payload.solverContractVersion !== GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION) {
      return {
        status: "engine-mismatch",
        label: "MTF sampling: provisional",
        solverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
        diagnostics: {
          warning: "Worker received a different geometric MTF solver contract.",
          expectedSolverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
          payloadSolverContractVersion: payload.solverContractVersion
        },
        comparisons: [],
        elapsedMs: Date.now() - started
      };
    }
    const surfaces = (payload.surfaceModel || payload.surfaces || [])
      .map((surface) => ({
        ...surface,
        x: toNumber(surface.x),
        radius: toNumber(surface.radius) || 0,
        semiDiameter: Math.max(0.001, toNumber(surface.semiDiameter) || 1),
        nBefore: toNumber(surface.nBefore) || 1,
        nAfter: toNumber(surface.nAfter) || 1,
        decenterY: toNumber(surface.decenterY) || 0,
        decenterZ: toNumber(surface.decenterZ) || 0,
        tiltY: toNumber(surface.tiltY) || 0,
        tiltZ: toNumber(surface.tiltZ) || 0,
        asphere: surface.asphere ? {
          ...surface.asphere,
          active: surface.asphere.active === true || asphereHasTerms(surface.asphere),
          enabled: surface.asphere.enabled === true || surface.asphere.active === true
        } : null,
        isStop: surface.isStop === true
      }))
      .filter((surface) => Number.isFinite(surface.x))
      .sort((a, b) => b.x - a.x);
    const modelSignature = surfaceSignature(surfaces);
    const features = surfaceFeatureFlags(surfaces);
    if (payload.expectedSurfaceSignature && payload.expectedSurfaceSignature !== modelSignature) {
      return {
        status: "engine-mismatch",
        label: "MTF sampling: provisional",
        solverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
        coreVersion: GEOMETRIC_MTF_CORE_VERSION,
        surfaceSignature: modelSignature,
        diagnostics: {
          warning: "Worker surface signature does not match the active optical model.",
          expectedSurfaceSignature: payload.expectedSurfaceSignature,
          surfaceSignature: modelSignature
        },
        comparisons: [],
        elapsedMs: Date.now() - started
      };
    }
    if (features.unsupportedByWorker) {
      return {
        status: "unsupported-geometry",
        label: "MTF sampling: provisional — worker model does not yet support active asphere/tilt geometry.",
        solverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
        coreVersion: GEOMETRIC_MTF_CORE_VERSION,
        surfaceSignature: modelSignature,
        diagnostics: {
          unsupportedReason: `Worker model does not yet support ${features.unsupportedFeatures.join(" and ")} geometry.`,
          unsupportedFeatures: features.unsupportedFeatures,
          surfaceFeatureFlags: features,
          warning: "Convergence was not calculated because the worker would otherwise solve a simplified optical model."
        },
        comparisons: [],
        fieldTargets: payload.fieldTargets || [],
        elapsedMs: Date.now() - started
      };
    }
    if (!surfaces.length) {
      return {
        status: "invalid",
        solverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
        coreVersion: GEOMETRIC_MTF_CORE_VERSION,
        surfaceSignature: modelSignature,
        comparisons: [],
        diagnostics: { warning: "No worker surface model supplied." },
        elapsedMs: Date.now() - started
      };
    }
    const apertureSurface = surfaces.find((surface) => surface.isStop) || surfaces[0];
    const rightMostSurfaceX = Math.max(...surfaces.map((surface) => surface.x));
    const setup = {
      orientation: payload.fieldOrientation || "tangential",
      apertureRadius: Math.max(0.001, toNumber(apertureSurface.semiDiameter) || (toNumber(payload.apertureDiameter) || 1) / 2),
      referenceX: apertureSurface.x,
      startX: rightMostSurfaceX + Math.max(20, Math.abs(toNumber(payload.totalTrackMm) || 0) * 0.15),
      imagePlaneX: Number.isFinite(toNumber(payload.imagePlaneX)) ? toNumber(payload.imagePlaneX) : 0,
      maxFieldAngleDegrees: payload.maxFieldAngleDegrees || 45
    };
    const frequencies = (payload.frequencies || [10, 30, 40, 50]).map(toNumber).filter(Number.isFinite);
    const maxFrequencyLpMm = Math.max(...frequencies, toNumber(payload.maxFrequencyLpMm) || 100);
    const fieldResults = (payload.fieldTargets || []).map((field) => convergenceForField(
      surfaces,
      setup,
      field,
      frequencies,
      maxFrequencyLpMm
    ));
    const comparisons = fieldResults.flatMap((result) => result.comparisons);
    const compared = comparisons.filter((item) => item.status === "compared");
    const interactiveHighDeltas = compared.map((item) => item.interactiveHighDelta).filter(Number.isFinite);
    const highReferenceDeltas = compared.map((item) => item.highReferenceDelta).filter(Number.isFinite);
    const maxInteractiveHighDelta = interactiveHighDeltas.length ? Math.max(...interactiveHighDeltas) : NaN;
    const maxHighReferenceDelta = highReferenceDeltas.length ? Math.max(...highReferenceDeltas) : NaN;
    const qualityDiagnostics = fieldResults.map((result) => ({
      fieldKey: result.field.key,
      fieldName: result.field.name,
      diagnostics: result.qualityDiagnostics
    }));
    const allDiagnostics = qualityDiagnostics.flatMap((field) => Object.values(field.diagnostics || {}));
    const clippedPresent = allDiagnostics.some((item) => item.clippedPupilWeight > 0.01);
    const minEffectiveReference = Math.min(...allDiagnostics.map((item) => item.effectiveSampleCount).filter(Number.isFinite));
    const allComparisonsCompleted = compared.length === (payload.fieldTargets || []).length * 2 * frequencies.length;
    const frequencyLimited = compared.some((item) => item.lsf?.kernel?.frequencyLimited);
    let status = "sampling-provisional";
    if (allDiagnostics.length && minEffectiveReference < 80) {
      status = "insufficient-clipped-pupil-coverage";
    } else if (
      allComparisonsCompleted
      && interactiveHighDeltas.every((delta) => delta <= 0.02)
      && highReferenceDeltas.every((delta) => delta <= 0.01)
      && !frequencyLimited
    ) {
      status = "converged";
    }
    const delta40 = Math.max(
      ...compared
        .filter((item) => item.frequency === 40)
        .flatMap((item) => [item.interactiveHighDelta, item.highReferenceDelta])
        .filter(Number.isFinite),
      NaN
    );
    return {
      status,
      label: status === "converged"
        ? "MTF sampling: converged"
        : status === "insufficient-clipped-pupil-coverage"
          ? "MTF sampling: insufficient clipped-pupil coverage"
          : "MTF sampling: provisional",
      fieldTargets: payload.fieldTargets || [],
      comparisons,
      maxInteractiveHighDelta,
      maxHighReferenceDelta,
      delta40,
      diagnostics: {
        solverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
        coreVersion: GEOMETRIC_MTF_CORE_VERSION,
        surfaceSignature: modelSignature,
        surfaceFeatureFlags: features,
        qualityDiagnostics,
        clippedPresent,
        minEffectiveSampleCount: Number.isFinite(minEffectiveReference) ? minEffectiveReference : NaN,
        allComparisonsCompleted,
        sharedLsfGridUsed: compared.every((item) => item.sharedLsfGridUsed),
        frequencyLimited,
        elapsedMs: Date.now() - started
      },
      solverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
      coreVersion: GEOMETRIC_MTF_CORE_VERSION,
      surfaceSignature: modelSignature,
      elapsedMs: Date.now() - started
    };
  };

  const normalizeSurfaceModel = (payload = {}) => (payload.surfaceModel || payload.surfaces || [])
    .map((surface) => ({
      ...surface,
      x: toNumber(surface.x),
      radius: toNumber(surface.radius) || 0,
      semiDiameter: Math.max(0.001, toNumber(surface.semiDiameter) || 1),
      nBefore: toNumber(surface.nBefore) || 1,
      nAfter: toNumber(surface.nAfter) || 1,
      decenterY: toNumber(surface.decenterY) || 0,
      decenterZ: toNumber(surface.decenterZ) || 0,
      tiltY: toNumber(surface.tiltY) || 0,
      tiltZ: toNumber(surface.tiltZ) || 0,
      asphere: surface.asphere ? {
        ...surface.asphere,
        active: surface.asphere.active === true || asphereHasTerms(surface.asphere),
        enabled: surface.asphere.enabled === true || surface.asphere.active === true
      } : null,
      isStop: surface.isStop === true
    }))
    .filter((surface) => Number.isFinite(surface.x))
    .sort((a, b) => b.x - a.x);

  const setupFromSurfaces = (surfaces, payload = {}) => {
    const apertureSurface = surfaces.find((surface) => surface.isStop) || surfaces[0];
    const rightMostSurfaceX = Math.max(...surfaces.map((surface) => surface.x));
    return {
      orientation: payload.fieldOrientation || "tangential",
      quality: mainPanelQuality(payload.quality),
      apertureRadius: Math.max(0.001, toNumber(apertureSurface?.semiDiameter) || (toNumber(payload.apertureDiameter) || 1) / 2),
      referenceX: apertureSurface?.x || 0,
      startX: rightMostSurfaceX + Math.max(20, Math.abs(toNumber(payload.totalTrackMm) || 0) * 0.15),
      imagePlaneX: Number.isFinite(toNumber(payload.imagePlaneX)) ? toNumber(payload.imagePlaneX) : 0,
      maxFieldAngleDegrees: payload.maxFieldAngleDegrees || 45
    };
  };

  const normalizeMainPanelPayload = (payload = {}) => {
    if (payload.solverContractVersion && payload.solverContractVersion !== GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION) {
      return {
        error: {
          status: "engine-mismatch",
          reason: "Worker received a different geometric MTF solver contract.",
          expectedSolverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
          payloadSolverContractVersion: payload.solverContractVersion
        }
      };
    }
    const surfaces = normalizeSurfaceModel(payload);
    const modelSignature = surfaceSignature(surfaces);
    if (payload.expectedSurfaceSignature && payload.expectedSurfaceSignature !== modelSignature) {
      return {
        surfaces,
        surfaceSignature: modelSignature,
        error: {
          status: "engine-mismatch",
          reason: "Worker surface signature does not match the active optical model.",
          expectedSurfaceSignature: payload.expectedSurfaceSignature,
          surfaceSignature: modelSignature
        }
      };
    }
    const features = surfaceFeatureFlags(surfaces);
    if (features.unsupportedByWorker) {
      return {
        surfaces,
        surfaceSignature: modelSignature,
        features,
        error: {
          status: "unsupported-geometry",
          reason: `Worker model does not yet support ${features.unsupportedFeatures.join(" and ")} geometry.`,
          unsupportedFeatures: features.unsupportedFeatures
        }
      };
    }
    if (!surfaces.length) {
      return {
        surfaces,
        surfaceSignature: modelSignature,
        features,
        error: {
          status: "invalid",
          reason: "No worker surface model supplied."
        }
      };
    }
    return { surfaces, surfaceSignature: modelSignature, features };
  };

  const calculateMainComparisons = (surfaces, setup, payload = {}) => {
    const requests = Array.isArray(payload.fieldRequests) ? payload.fieldRequests : [];
    return requests.map((request) => {
      const current = calculateMainGeometricLsfMtf(surfaces, setup, {
        ...payload,
        ...request,
        maxFrequencyLpMm: payload.maxFrequencyLpMm,
        frequencyStepLpMm: payload.frequencyStepLpMm || 5,
        quality: payload.quality,
        spectralLineKey: request.spectralLineKey || payload.wavelength?.spectralLineKey || "d",
        wavelengthNm: request.wavelengthNm || payload.wavelength?.wavelengthNm || 587.6
      });
      return {
        fieldKey: request.fieldKey,
        fieldName: request.fieldName || "Field",
        fieldAngleDegrees: request.fieldAngleDegrees,
        spectralLineKey: request.spectralLineKey || "d",
        current,
        best: null,
        bestFocusDeferred: true,
        enginesMatch: true,
        comparisonUnavailableReason: ""
      };
    });
  };

  const calculateManufacturerFieldData = (surfaces, setup, payload = {}, request = {}) => {
    const maxHeight = Math.max(0, toNumber(request.maxImageHeightMm) || 0);
    const sampleCount = Math.max(3, Math.round(toNumber(request.sampleCount) || 9));
    const frequencies = (request.frequencies || [10, 30]).map(toNumber).filter(Number.isFinite);
    const warnings = new Set();
    const sampleForTarget = (targetImageHeight, index) => {
      const solved = solveFieldAngleForImageHeight(surfaces, setup, targetImageHeight);
      const result = solved.status === "solved" || targetImageHeight <= 1e-6
        ? calculateMainGeometricLsfMtf(surfaces, setup, {
          ...payload,
          fieldAngleDegrees: solved.fieldAngleDegrees || 0,
          fieldName: `${Number(targetImageHeight).toFixed(1)} mm`,
          fieldKey: `height-${index}`,
          maxFrequencyLpMm: payload.maxFrequencyLpMm,
          frequencyStepLpMm: payload.frequencyStepLpMm || 5,
          quality: payload.quality
        })
        : null;
      if (solved.status !== "solved") warnings.add(`Target image height ${targetImageHeight.toFixed(2)} mm could not be solved (${solved.status}).`);
      const imageHeight = Number.isFinite(solved.imageHeight) ? solved.imageHeight : NaN;
      const isSolved = solved.status === "solved" && Number.isFinite(imageHeight) && result;
      return {
        targetImageHeight,
        imageHeight,
        fieldAngleDegrees: solved.fieldAngleDegrees,
        solvedFromChiefRay: solved.status === "solved",
        solveStatus: solved.status,
        residualMm: solved.residualMm,
        result: result || {
          status: solved.status,
          engine: "geometricLsfFft",
          engineLabel: "Geometric LSF/FFT preview — not physical wavefront MTF",
          quality: mainPanelQuality(payload.quality),
          tangential: { frequencies: [], values: [] },
          sagittal: { frequencies: [], values: [] }
        },
        readouts: Object.fromEntries(frequencies.map((frequency) => [
          frequency,
          {
            sagittal: isSolved ? mtfValueAtFrequency(result.sagittal, frequency) : NaN,
            tangential: isSolved ? mtfValueAtFrequency(result.tangential, frequency) : NaN
          }
        ]))
      };
    };
    const targets = Array.from({ length: sampleCount }, (_, index) => (maxHeight * index) / (sampleCount - 1));
    const samples = targets.map(sampleForTarget).sort((left, right) => (
      (Number.isFinite(left.imageHeight) ? left.imageHeight : Number.POSITIVE_INFINITY)
      - (Number.isFinite(right.imageHeight) ? right.imageHeight : Number.POSITIVE_INFINITY)
    ));
    const solvedEdgeMm = samples.reduce((edge, sample) => (
      sample.solveStatus === "solved" && Number.isFinite(sample.imageHeight)
        ? Math.max(edge, sample.imageHeight)
        : edge
    ), 0);
    return {
      maxImageHeightMm: maxHeight,
      selectedSensorCornerMm: maxHeight,
      solvedEdgeMm,
      reachedSensorCorner: samples.some((sample) => sample.solveStatus === "solved" && Number.isFinite(sample.imageHeight) && Math.abs(sample.imageHeight - maxHeight) <= 0.02),
      sensor: request.sensor || null,
      frequencies,
      engine: "geometricLsfFft",
      engineLabel: "Geometric LSF/FFT preview — not physical wavefront MTF",
      quality: mainPanelQuality(payload.quality),
      focusPolicy: request.focusPolicy || "fixed",
      samples,
      source: "worker · iterated chief-ray intercept · geometric LSF/FFT",
      imagePlaneX: setup.imagePlaneX,
      fieldSampleCount: samples.length,
      warnings: [...warnings]
    };
  };

  const calculateApertureSweep = (payload = {}, baseSetup, baseSurfaces) => {
    const request = payload.apertureSweepRequest || {};
    const options = Array.isArray(request.options) ? request.options : [];
    const field = request.field || { key: "center", name: "Centre", angle: 0 };
    return options.map((option) => {
      const surfaces = option.surfaceModel ? normalizeSurfaceModel({ surfaceModel: option.surfaceModel }) : baseSurfaces;
      const setup = setupFromSurfaces(surfaces, {
        ...payload,
        imagePlaneX: Number.isFinite(toNumber(option.imagePlaneX)) ? toNumber(option.imagePlaneX) : baseSetup.imagePlaneX
      });
      const result = calculateMainGeometricLsfMtf(surfaces, setup, {
        ...payload,
        fieldKey: field.key,
        fieldName: field.name,
        fieldAngleDegrees: field.angle || 0,
        maxFrequencyLpMm: payload.maxFrequencyLpMm,
        frequencyStepLpMm: payload.frequencyStepLpMm || 5,
        quality: payload.quality
      });
      const fieldData = calculateManufacturerFieldData(surfaces, setup, payload, {
        ...(payload.manufacturerRequest || {}),
        ...(option.manufacturerRequest || {}),
        sampleCount: request.chartMode === "field" ? 7 : (payload.manufacturerRequest?.sampleCount || 9),
        focusPolicy: option.focusPolicy || "fixed"
      });
      return {
        key: option.key,
        label: option.label,
        fNumber: option.fNumber ?? null,
        requestedFNumber: option.requestedFNumber ?? option.fNumber ?? null,
        physicalStopDiameter: option.physicalStopDiameter,
        imagePlaneX: setup.imagePlaneX,
        bestFocusStable: option.bestFocusStable,
        refocusAtSearchLimit: option.refocusAtSearchLimit === true,
        focusPolicy: option.focusPolicy || "fixed",
        fieldData,
        result,
        retraced: true,
        cachedAtRevision: payload.cacheRevision || 0,
        source: "worker"
      };
    });
  };

  const calculateGeometricLsfMainPanel = (payload = {}) => {
    const started = Date.now();
    const normalized = normalizeMainPanelPayload(payload);
    if (normalized.error) {
      return {
        status: normalized.error.status,
        solverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
        surfaceSignature: normalized.surfaceSignature || "",
        warnings: [normalized.error.reason],
        diagnostics: normalized.error,
        source: "worker",
        comparisons: [],
        activeResults: [],
        manufacturerFieldData: null,
        apertureSweepResults: [],
        elapsedMs: Date.now() - started
      };
    }
    const { surfaces, surfaceSignature: modelSignature, features } = normalized;
    const setup = setupFromSurfaces(surfaces, payload);
    const comparisons = calculateMainComparisons(surfaces, setup, payload);
    const activeResults = comparisons.map((comparison) => comparison.current);
    const manufacturerFieldData = payload.manufacturerRequest
      ? calculateManufacturerFieldData(surfaces, setup, payload, payload.manufacturerRequest)
      : null;
    const apertureSweepResults = calculateApertureSweep(payload, setup, surfaces);
    const warnings = new Set();
    activeResults.forEach((result) => (result.warnings || []).forEach((warning) => warnings.add(warning)));
    (manufacturerFieldData?.warnings || []).forEach((warning) => warnings.add(warning));
    apertureSweepResults.forEach((item) => {
      (item.result?.warnings || []).forEach((warning) => warnings.add(warning));
      (item.fieldData?.warnings || []).forEach((warning) => warnings.add(warning));
    });
    return {
      status: "valid",
      source: "worker",
      solverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
      coreVersion: GEOMETRIC_MTF_CORE_VERSION,
      workerVersion: payload.workerVersion || "",
      surfaceSignature: modelSignature,
      surfaceFeatureFlags: features,
      comparisons,
      activeResults,
      manufacturerFieldData,
      apertureSweepResults,
      warnings: [...warnings],
      rayCount: payload.rayCount,
      maxFrequencyLpMm: payload.maxFrequencyLpMm || 100,
      diagnostics: {
        elapsedMs: Date.now() - started,
        surfaceSignature: modelSignature,
        solverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
        coreVersion: GEOMETRIC_MTF_CORE_VERSION,
        workerEligible: true
      },
      elapsedMs: Date.now() - started
    };
  };

  const api = {
    GEOMETRIC_MTF_CORE_VERSION,
    GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
    surfaceFeatureFlags,
    surfaceSignature,
    calculateGeometricLsfMainPanel,
    calculateGeometricLsfConvergence,
    QUALITY_PROFILES
  };

  if (typeof self !== "undefined") self.geometricMtfCore = api;
  if (typeof globalThis !== "undefined") globalThis.geometricMtfCore = api;
})();

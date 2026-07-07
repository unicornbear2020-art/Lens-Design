/* eslint-env worker */

const GEOMETRIC_MTF_WORKER_VERSION = "20260707-ipad-system-result-1";
const GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION = "geometric-lsf-contract-20260630-1";

try {
  importScripts("geometric-mtf-core.js?v=20260707-ipad-system-result-1");
} catch (error) {
  self.geometricMtfCoreLoadError = error;
}

self.onmessage = (event) => {
  const payload = event.data || {};
  const { requestId, task } = payload;

  if (!["geometric-lsf-convergence", "geometric-lsf-main-panel"].includes(task)) {
    self.postMessage({
      requestId,
      status: "error",
      workerVersion: GEOMETRIC_MTF_WORKER_VERSION,
      coreVersion: self.geometricMtfCore?.GEOMETRIC_MTF_CORE_VERSION || "",
      solverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
      error: `Unsupported geometric MTF worker task: ${task || "unknown"}`
    });
    return;
  }

  const coreFunction = task === "geometric-lsf-main-panel"
    ? self.geometricMtfCore?.calculateGeometricLsfMainPanel
    : self.geometricMtfCore?.calculateGeometricLsfConvergence;

  if (!coreFunction) {
    self.postMessage({
      requestId,
      status: "error",
      workerVersion: GEOMETRIC_MTF_WORKER_VERSION,
      coreVersion: self.geometricMtfCore?.GEOMETRIC_MTF_CORE_VERSION || "",
      error: self.geometricMtfCoreLoadError?.message || "geometric-mtf-core.js did not load."
    });
    return;
  }

  try {
    const result = coreFunction(payload);
    self.postMessage({
      requestId,
      status: "complete",
      workerVersion: GEOMETRIC_MTF_WORKER_VERSION,
      coreVersion: self.geometricMtfCore?.GEOMETRIC_MTF_CORE_VERSION || "",
      solverContractVersion: GEOMETRIC_MTF_SOLVER_CONTRACT_VERSION,
      surfaceSignature: result?.surfaceSignature || payload.expectedSurfaceSignature || "",
      result
    });
  } catch (error) {
    self.postMessage({
      requestId,
      status: "error",
      workerVersion: GEOMETRIC_MTF_WORKER_VERSION,
      coreVersion: self.geometricMtfCore?.GEOMETRIC_MTF_CORE_VERSION || "",
      error: error?.message || String(error)
    });
  }
};

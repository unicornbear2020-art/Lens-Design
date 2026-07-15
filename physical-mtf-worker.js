/* eslint-env worker */

const PHYSICAL_MTF_WORKER_VERSION = "20260715-physical-mtf-focus-response-1";

try {
  importScripts("physical-mtf-core.js?v=20260715-physical-mtf-focus-response-1");
} catch (error) {
  self.physicalMtfCoreLoadError = error;
}

self.onmessage = (event) => {
  const payload = event.data || {};
  const { requestId, gridSize, task } = payload;
  const core = self.physicalMtfCore;
  try {
    if (!core) {
      throw self.physicalMtfCoreLoadError || new Error("Physical MTF core is unavailable.");
    }
    const result = task === "lens-complex-pupil"
      ? core.calculateLensComplexPupilMtf(payload)
      : task === "canonical-aberration-validation"
        ? core.calculateCanonicalAberrationValidation()
        : core.calculateIdealCircularPupilValidation(gridSize);
    self.postMessage({
      requestId,
      status: result.status === "invalid" || result.status === "prototype-failed" ? "error" : "complete",
      workerVersion: PHYSICAL_MTF_WORKER_VERSION,
      coreVersion: core.PHYSICAL_MTF_CORE_VERSION,
      task: task || "ideal-circular-pupil",
      result
    });
  } catch (error) {
    self.postMessage({
      requestId,
      status: "error",
      workerVersion: PHYSICAL_MTF_WORKER_VERSION,
      coreVersion: core?.PHYSICAL_MTF_CORE_VERSION || "",
      task: task || "ideal-circular-pupil",
      result: { status: "invalid", warning: error?.message || "Worker calculation failed." }
    });
  }
};

# Aperture Stop Audit

This audit documents how each patent-based preset resolves its aperture stop in the app. It does not claim production-lens certainty unless the local prescription metadata explicitly identifies a stop surface. When a source is listed as estimated, the app uses it only as an editable optical-analysis default.

## Resolution Rules

1. Manual user settings have priority: surface number, distance from sensor, or distance from the first/front surface.
2. Surface prescriptions in Auto mode use `apertureStopSpec` metadata first.
3. If no `apertureStopSpec` exists, legacy `surface.isStop` rows are used.
4. If neither exists, the app estimates the central non-cemented optical air gap and labels the result as unverified.
5. The stop location affects ray tracing, f-number, illumination, spot, MTF, wavefront, and coverage calculations. It does not resize lens elements.

## Patent Preset Stop Sources

| Preset key | Display name | Resolved stop type | Surface / air gap | Source level | Confidence | Source URL | Reason | Production match confirmed? |
|---|---|---:|---|---|---|---|---|---|
| `zunow50F11Us2715354` | Zunow 50mm f/1.1 — US2715354 | Auto central air-gap estimate | Computed at runtime | Estimated | Unverified | https://patents.google.com/patent/US2715354A/en | No stop plane is confirmed in the current local transcription. | No |
| `gaussF2Us4123144Ex1` | Patent Gauss f/2 — US4123144 Example 1 | Patent stop surface | Surface 6 | Patent | Verified | https://patents.google.com/patent/US4123144A/en | The entered prescription table marks surface 6 as STOP. | Patent example only |
| `leicaSummilux50F14Us3291553` | Leica Summilux 50mm f/1.4 — US3291553A | Auto central air-gap estimate | Computed at runtime | Estimated | Unverified | https://patents.google.com/patent/US3291553A/en | No stop plane is confirmed in the current local transcription. | No |
| `leicaSummicronC40F2De2222892Ex3` | Leica Summicron-C 40mm f/2 — DE2222892 Example 3 | Auto central air-gap estimate | Computed at runtime | Estimated | Unverified | https://patents.google.com/patent/DE2222892A1/en | No stop plane is confirmed in the current local transcription. | No |
| `canonDreamLens50F095JpSho39_10178` | Canon 50mm f/0.95 "Dream Lens" — Tokkosho 39-10178 | Estimated patent-layout stop | Gap S7-S8 midpoint | Estimated | Unverified |  | Patent gives f/0.95 but not a numeric stop coordinate; the central S3 region is used as an editable estimate. | No |
| `zeissSonnar50F15Us2186621Ex1` | Zeiss Sonnar 50mm f/1.5 — US2186621 Example 1 | Auto central air-gap estimate | Computed at runtime | Estimated | Unverified | https://patents.google.com/patent/US2186621A/en | No stop plane is confirmed in the current local transcription. | No |
| `zeissSonnarType50F14Us2600610Ex3` | Zeiss Sonnar-type 50mm f/1.4 — US2600610 Example III | Estimated patent-layout stop | Gap S6-S7 midpoint | Estimated | Probable | https://patents.google.com/patent/US2600610A/en | The layout places the stop in the air gap between the front and rear components, but the exact fraction is not stated in the local numeric prescription. | No |
| `zeissSonnar40F28Us3994576Ex1` | Zeiss Sonnar 40mm f/2.8 — US3994576 Example 1 | Auto central air-gap estimate | Computed at runtime | Estimated | Unverified | https://patents.google.com/patent/US3994576A/en | No stop plane is confirmed in the current local transcription. | No |
| `zeissBiotar50F14Us1786916Ex2` | Carl Zeiss Biotar 50mm f/1.4 — US1786916 Example 2 | Estimated patent-layout stop | Gap S5-S6 midpoint | Estimated | Probable | https://patents.google.com/patent/US1786916A/en | The central air space between the two cemented halves is used as an editable estimate. | No |
| `sonyZeissPlanar50F14Us20140071331Ex4` | Sony Zeiss Planar 50mm f/1.4 ZA — US20140071331 Example 4 | Patent stop surface | Surface 8 | Patent | Verified | https://patents.google.com/patent/US20140071331A1/en | The entered prescription table includes stop surface 8. | Patent example only |
| `sonySonnarFe55F18Us20150092100Ex1` | Sony Sonnar T* FE 55mm f/1.8 ZA — US20150092100 Example 1 | Patent stop surface | Surface 7 | Patent | Verified | https://patents.google.com/patent/US20150092100A1/en | The entered prescription table includes stop surface 7. | Patent example only |
| `sonyFe24F14Wo2019073744Ex1` | Sony FE 24mm f/1.4 GM — WO2019-073744 Example 1 | Auto central air-gap estimate | Computed at runtime | Estimated | Unverified | https://patents.google.com/patent/WO2019073744A1/en | No stop plane is confirmed in the current local transcription. | No |
| `sonyFe70200F4Us20150226945Ex1` | Sony FE 70-200mm f/4 G OSS — US20150226945 Example 1 | Auto central air-gap estimate | Computed at runtime | Estimated | Unverified | https://patents.google.com/patent/US20150226945A1/en | No stop plane is confirmed in the current local transcription. | No |
| `nikonLargeApertureWideAngleUs3576360` | Nikon large aperture wide angle — US3576360A | Auto central air-gap estimate | Computed at runtime | Estimated | Unverified | https://patents.google.com/patent/US3576360A/en | No stop plane is confirmed in the current local transcription. | No |
| `nikonF12Us3738736` | Nikon f/1.2 lens — US3738736A | Auto central air-gap estimate | Computed at runtime | Estimated | Unverified | https://patents.google.com/patent/US3738736A/en | No stop plane is confirmed in the current local transcription. | No |
| `olympusF12Us4099843` | Olympus f/1.2 lens — US4099843A | Auto central air-gap estimate | Computed at runtime | Estimated | Unverified | https://patents.google.com/patent/US4099843A/en | No stop plane is confirmed in the current local transcription. | No |

## Known Limitations

- Several local patent transcriptions do not include an explicit aperture stop row. Those presets are intentionally marked estimated instead of silently using the front of the lens group.
- Estimated stops are analysis defaults only. They remain editable through the Aperture & Field Setup panel.
- Production match is not confirmed by this audit unless separately supported by production service data.

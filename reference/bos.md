# Bronchiolitis obliterans syndrome after lung transplants

A dataset containing histories of bronchiolitis obliterans syndrome
(BOS) from lung transplant recipients. BOS is a chronic decline in lung
function, often observed after lung transplantation.

## Format

A data frame containing a sequence of observed or censored transitions
to the next stage of severity or death. It is grouped by patient and
includes histories of 204 patients. All patients start in state 1 (no
BOS) at six months after transplant, and may subsequently develop BOS or
die.

`bosms3` contains the data for a three-state model: no BOS, BOS or
death. `bosms4` uses a four-state representation: no BOS, mild BOS,
moderate/severe BOS or death.

|          |           |                                                                                                                                                                           |
|----------|-----------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `id`     | (numeric) | Patient identification number                                                                                                                                             |
| `from`   | (numeric) | Observed starting state of the transition                                                                                                                                 |
| `to`     | (numeric) | Observed or potential ending state of the transition                                                                                                                      |
| `Tstart` | (numeric) | Time at the start of the interval                                                                                                                                         |
| `Tstop`  | (numeric) | Time at the end of the interval                                                                                                                                           |
| `time`   | (numeric) | Time difference `Tstart`-`Tstop`                                                                                                                                          |
| `status` | (numeric) | 1 if the transition to state `to` was observed, or 0 if the transition to state `to` was censored (for example, if the patient was observed to move to a competing state) |
| `trans`  | (factor)  | Number of the transition `from`-`to` in the set of all `ntrans` allowed transitions, numbered from 1 to `ntrans`.                                                         |

## Source

Papworth Hospital, U.K.

## Details

The entry time of each patient into each stage of BOS was estimated by
clinicians, based on their history of lung function measurements and
acute rejection and infection episodes. BOS is only assumed to occur
beyond six months after transplant. In the first six months the function
of each patient's new lung stabilises. Subsequently BOS is diagnosed by
comparing the lung function against the "baseline" value.

The same data are provided in the msm package, but in the native format
of msm to allow Markov models to be fitted. In flexsurv, much more
flexible models can be fitted.

## References

Heng. D. et al. (1998). Bronchiolitis Obliterans Syndrome: Incidence,
Natural History, Prognosis, and Risk Factors. Journal of Heart and Lung
Transplantation 17(12)1255–1263.

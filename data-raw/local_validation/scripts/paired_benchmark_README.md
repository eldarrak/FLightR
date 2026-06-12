# Paired Local Benchmark Runner

This local-only runner compares:

- baseline checkout: `D:/GitHub/FLightR_fixed_baseline`
- s2/current checkout: `D:/GitHub/FLightR`

Private data stay in `D:/GitHub/FLightR/data-raw/local_validation/private_data/`.
All benchmark outputs and logs are written under
`D:/GitHub/FLightR/data-raw/local_validation/outputs/`.

## Dry Run From R

From the FLightR checkout:

```r
source("data-raw/local_validation/scripts/run_paired_benchmark_once.R")
```

The R runner is dry-run by default. It prints the labels and settings it would run.

To run one selected version from R, edit the settings near the top of
`run_paired_benchmark_once.R`, or set environment variables such as:

```r
Sys.setenv(
  FLIGHTR_PAIR_ACTUALLY_RUN = "true",
  FLIGHTR_PAIR_RUN_BASELINE = "true",
  FLIGHTR_PAIR_RUN_S2 = "false"
)
source("data-raw/local_validation/scripts/run_paired_benchmark_once.R")
```

## Actual Parallel Run From PowerShell

Edit settings near the top of:

```text
D:/GitHub/FLightR/data-raw/local_validation/scripts/run_paired_benchmark_parallel.ps1
```

Defaults:

- `nParticles = 100000`
- `prerun threads per run = 4, PF threads per run = 1`
- `seed = 123`
- `pair_id = "pair1"`

Run:

```powershell
powershell -ExecutionPolicy Bypass -File D:/GitHub/FLightR/data-raw/local_validation/scripts/run_paired_benchmark_parallel.ps1
```

The PowerShell runner starts baseline and s2 as separate `Rscript` processes close
together, waits for both to finish, and writes logs to:

```text
D:/GitHub/FLightR/data-raw/local_validation/outputs/logs/
```

## Summarize Results

From R:

```r
source("data-raw/local_validation/scripts/summarize_paired_benchmark.R")
```

This writes:

- `paired_runtime_summary.csv`
- `paired_runtime_comparison.csv`

Both files are saved under `data-raw/local_validation/outputs/`.

## Notes

Use `prerun threads per run = 4, PF threads per run = 1` by default because running two 8-thread jobs at the
same time can oversubscribe the machine.

Do not use `1e6` by default. First confirm the paired runner works at `1e5`;
then manually increase `nParticles` if the runtime is acceptable.

Sleep interruptions make absolute runtime unreliable. Paired parallel runs make
sleep time more similar between versions, but the two jobs still compete for CPU,
memory bandwidth, and thermal limits. If speed differences are small, follow up
with a no-sleep sequential benchmark.

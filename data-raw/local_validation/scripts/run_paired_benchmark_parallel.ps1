# Local-only paired parallel benchmark runner for FLightR.
# Starts baseline and s2 runs at nearly the same time.
# All logs and generated runners are written under data-raw/local_validation/outputs/.

$ErrorActionPreference = "Stop"

# ----------------------------
# Editable settings
# ----------------------------

$nParticles = 100000
$threads_per_run = 4
$seed = 123
$pair_id = "pair1"

# Leave empty to auto-detect.
$RscriptPath = ""

$BaselineRepo = "D:/GitHub/FLightR_fixed_baseline"
$S2Repo = "D:/GitHub/FLightR"
$ScriptDir = "D:/GitHub/FLightR/data-raw/local_validation/scripts"
$OutputDir = "D:/GitHub/FLightR/data-raw/local_validation/outputs"

# ----------------------------
# Helpers
# ----------------------------

function Get-ParticleLabel([int]$n) {
    $pow = [Math]::Log10($n)
    if ([Math]::Abs($pow - [Math]::Round($pow)) -lt 1e-10) {
        return "1e$([int][Math]::Round($pow))"
    }
    return "$n"
}

function Find-Rscript {
    param([string]$ManualPath)
    if ($ManualPath -and (Test-Path -LiteralPath $ManualPath)) {
        return $ManualPath
    }
    $candidates = @(
        "C:/Program Files/R/R-4.6.0/bin/Rscript.exe",
        "C:/Program Files/R/R-4.5.0/bin/Rscript.exe",
        "C:/Program Files/R/R-4.4.0/bin/Rscript.exe",
        "C:/Program Files/R/R-4.3.0/bin/Rscript.exe"
    )
    foreach ($candidate in $candidates) {
        if (Test-Path -LiteralPath $candidate) {
            return $candidate
        }
    }
    $cmd = Get-Command Rscript.exe -ErrorAction SilentlyContinue
    if ($cmd) {
        return $cmd.Source
    }
    throw "Could not find Rscript.exe. Set `$RscriptPath near the top of this script."
}

function New-RunnerFile {
    param(
        [string]$Version,
        [string]$RepoPath,
        [string]$RunLabel,
        [string]$RunnerPath,
        [int]$ParallelPort
    )

    $runOnce = Join-Path $ScriptDir "run_paired_benchmark_once.R"
    $runBaseline = ([string]($Version -eq "baseline")).ToLowerInvariant()
    $runS2 = ([string]($Version -eq "s2")).ToLowerInvariant()
    $content = @"
Sys.setenv(
  FLIGHTR_PAIR_NPARTICLES = "$nParticles",
  FLIGHTR_PAIR_THREADS_PER_RUN = "$threads_per_run",
  FLIGHTR_PAIR_SEED = "$seed",
  FLIGHTR_PAIR_ID = "$pair_id",
  FLIGHTR_PAIR_ACTUALLY_RUN = "true",
  FLIGHTR_PAIR_RERUN_EXISTING = "false",
  FLIGHTR_PAIR_RUN_BASELINE = "$runBaseline",
  FLIGHTR_PAIR_RUN_S2 = "$runS2",
  R_PARALLEL_PORT = "$ParallelPort"
)
setwd("$($RepoPath -replace "\\", "/")")
source("$($runOnce -replace "\\", "/")")
"@
    $utf8NoBom = New-Object System.Text.UTF8Encoding($false)
    [System.IO.File]::WriteAllText($RunnerPath, $content, $utf8NoBom)
}

function Start-PairedRun {
    param(
        [string]$Version,
        [string]$RepoPath,
        [string]$RunLabel,
        [string]$Rscript,
        [int]$ParallelPort
    )

    if (!(Test-Path -LiteralPath $RepoPath)) {
        throw "Repository path not found for ${Version}: $RepoPath"
    }

    $tempDir = Join-Path $OutputDir "temp_runners"
    $logDir = Join-Path $OutputDir "logs"
    New-Item -ItemType Directory -Path $tempDir -Force | Out-Null
    New-Item -ItemType Directory -Path $logDir -Force | Out-Null

    $runnerPath = Join-Path $tempDir "${Version}_${RunLabel}.R"
    $stdoutPath = Join-Path $logDir "${Version}_${RunLabel}.log"
    $stderrPath = Join-Path $logDir "${Version}_${RunLabel}.err.log"

    New-RunnerFile -Version $Version -RepoPath $RepoPath -RunLabel $RunLabel -RunnerPath $runnerPath -ParallelPort $ParallelPort

    Write-Host "Starting ${Version}: $RunLabel"
    Write-Host "  R_PARALLEL_PORT: $ParallelPort"
    Write-Host "  log: $stdoutPath"
    Write-Host "  err: $stderrPath"

    $process = Start-Process `
        -FilePath $Rscript `
        -ArgumentList @($runnerPath) `
        -WorkingDirectory $RepoPath `
        -RedirectStandardOutput $stdoutPath `
        -RedirectStandardError $stderrPath `
        -WindowStyle Hidden `
        -PassThru

    [PSCustomObject]@{
        Version = $Version
        RunLabel = $RunLabel
        Process = $process
        Stdout = $stdoutPath
        Stderr = $stderrPath
        Runner = $runnerPath
    }
}

# ----------------------------
# Run
# ----------------------------

New-Item -ItemType Directory -Path $OutputDir -Force | Out-Null

$Rscript = Find-Rscript -ManualPath $RscriptPath
$particleLabel = Get-ParticleLabel -n $nParticles
$baselineLabel = "baseline_halfyear_${particleLabel}_${pair_id}_seed${seed}"
$s2Label = "s2_halfyear_${particleLabel}_${pair_id}_seed${seed}"

Write-Host "Rscript: $Rscript"
Write-Host "nParticles=$nParticles threads_per_run=$threads_per_run seed=$seed pair_id=$pair_id"
Write-Host "Outputs: $OutputDir"

$runs = @()
$runs += Start-PairedRun -Version "baseline" -RepoPath $BaselineRepo -RunLabel $baselineLabel -Rscript $Rscript -ParallelPort 12101
$runs += Start-PairedRun -Version "s2" -RepoPath $S2Repo -RunLabel $s2Label -Rscript $Rscript -ParallelPort 12102

foreach ($run in $runs) {
    $run.Process.WaitForExit()
}

Write-Host ""
Write-Host "Paired benchmark finished."
foreach ($run in $runs) {
    Write-Host "$($run.Version) exit=$($run.Process.ExitCode) label=$($run.RunLabel)"
    Write-Host "  log: $($run.Stdout)"
    Write-Host "  err: $($run.Stderr)"
}

param(
  [string]$Source = "src/sec4/kelbm.c",
  [string]$OutputDir = "outputs/sec4/kelbm",
  [switch]$SkipPlot
)

$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir

Push-Location $repoRoot
try {
  & (Join-Path $scriptDir "build_one.cmd") $Source
  if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
  }

  $exePath = Join-Path $repoRoot "build/bin/kelbm.exe"
  if (-not (Test-Path $exePath)) {
    throw "Executable was not generated: $exePath"
  }

  $runDir = Join-Path $repoRoot $OutputDir
  New-Item -ItemType Directory -Path $runDir -Force | Out-Null

  Write-Host "Running kelbm.exe in $runDir"
  Push-Location $runDir
  try {
    & $exePath
    if ($LASTEXITCODE -ne 0) {
      throw "kelbm.exe exited with code $LASTEXITCODE"
    }
  }
  finally {
    Pop-Location
  }

  if (-not $SkipPlot) {
    $venvPython = Join-Path $repoRoot ".venv/Scripts/python.exe"
    $python = if (Test-Path $venvPython) { $venvPython } else { "python" }
    Write-Host "Regenerating plots with $python"
    & $python (Join-Path $scriptDir "plot_kelbm_compare_theory.py")
    & $python (Join-Path $scriptDir "plot_kelbm_contour.py")
    & $python (Join-Path $scriptDir "plot_kelbm_centerline.py")
    & $python (Join-Path $scriptDir "plot_kelbm_output.py")
  }
}
finally {
  Pop-Location
}

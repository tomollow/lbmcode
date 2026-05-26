param(
  [switch]$SkipBuild,
  [switch]$SkipRun,
  [switch]$SkipPlot
)

$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir

Push-Location $repoRoot
try {
  $exePath = Join-Path $repoRoot "build/bin/lbmcm_karman.exe"
  $runDir  = Join-Path $repoRoot "outputs/sec4/lbmcm_karman"

  if (-not $SkipBuild) {
    Write-Host "Building src/sec4/lbmcm_karman.c"
    & (Join-Path $scriptDir "build_one.cmd") "src/sec4/lbmcm_karman.c"
    if ($LASTEXITCODE -ne 0) { throw "build_one.cmd exited with code $LASTEXITCODE" }
    if (-not (Test-Path $exePath)) { throw "Executable was not generated: $exePath" }
  }

  if (-not $SkipRun) {
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    Write-Host "Running lbmcm_karman.exe in $runDir"
    Push-Location $runDir
    try {
      & $exePath
      if ($LASTEXITCODE -ne 0) { throw "lbmcm_karman.exe exited with code $LASTEXITCODE" }
    } finally { Pop-Location }
  }

  if (-not $SkipPlot) {
    $venvPython = Join-Path $repoRoot ".venv/Scripts/python.exe"
    $python = if (Test-Path $venvPython) { $venvPython } else { "python" }
    Write-Host "Regenerating plot with $python"
    & $python (Join-Path $scriptDir "plot_lbmcm_karman_spectrum.py")
    if ($LASTEXITCODE -ne 0) { throw "plot script exited with code $LASTEXITCODE" }
  }
} finally { Pop-Location }

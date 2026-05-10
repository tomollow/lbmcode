param(
  [switch]$SkipPlot,
  [switch]$PureOnly,
  [switch]$KepsOnly
)

$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir

$variants = @()
if (-not $KepsOnly) {
  $variants += @{ Name = "taylor_green";      Source = "src/sec4/taylor_green.c";      Exe = "taylor_green.exe" }
}
if (-not $PureOnly) {
  $variants += @{ Name = "taylor_green_keps"; Source = "src/sec4/taylor_green_keps.c"; Exe = "taylor_green_keps.exe" }
}

Push-Location $repoRoot
try {
  foreach ($v in $variants) {
    & (Join-Path $scriptDir "build_one.cmd") $v.Source
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    $exePath = Join-Path $repoRoot "build/bin/$($v.Exe)"
    if (-not (Test-Path $exePath)) {
      throw "Executable was not generated: $exePath"
    }

    $runDir = Join-Path $repoRoot "outputs/sec4/$($v.Name)"
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null

    Write-Host "Running $($v.Exe) in $runDir"
    Push-Location $runDir
    try {
      & $exePath
      if ($LASTEXITCODE -ne 0) {
        throw "$($v.Exe) exited with code $LASTEXITCODE"
      }
    }
    finally { Pop-Location }
  }

  if (-not $SkipPlot) {
    $venvPython = Join-Path $repoRoot ".venv/Scripts/python.exe"
    $python = if (Test-Path $venvPython) { $venvPython } else { "python" }
    Write-Host "Regenerating plots with $python"
    if (-not $KepsOnly -and -not $PureOnly) {
      # Decay plot needs both runs
      & $python (Join-Path $scriptDir "plot_taylor_green_decay.py")
    }
    if (-not $KepsOnly) {
      & $python (Join-Path $scriptDir "plot_taylor_green_snapshots.py") "pure"
    }
    if (-not $PureOnly) {
      & $python (Join-Path $scriptDir "plot_taylor_green_snapshots.py") "keps"
    }
  }
}
finally { Pop-Location }

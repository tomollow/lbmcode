param(
  [switch]$SkipPlot,
  [switch]$PureOnly,
  [switch]$KepsOnly,
  [switch]$LesOnly
)

$ErrorActionPreference = "Stop"

$onlyCount = @($PureOnly, $KepsOnly, $LesOnly | Where-Object { $_ }).Count
if ($onlyCount -gt 1) {
  throw "Specify at most one of -PureOnly / -KepsOnly / -LesOnly"
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir

$variants = @()
if (-not $KepsOnly -and -not $LesOnly) {
  $variants += @{ Name = "kelvin_helmholtz";      Source = "src/sec4/kelvin_helmholtz.c";      Exe = "kelvin_helmholtz.exe" }
}
if (-not $PureOnly -and -not $LesOnly) {
  $variants += @{ Name = "kelvin_helmholtz_keps"; Source = "src/sec4/kelvin_helmholtz_keps.c"; Exe = "kelvin_helmholtz_keps.exe" }
}
if (-not $PureOnly -and -not $KepsOnly) {
  $variants += @{ Name = "kelvin_helmholtz_les";  Source = "src/sec4/kelvin_helmholtz_les.c";  Exe = "kelvin_helmholtz_les.exe" }
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
    if ($onlyCount -eq 0) {
      & $python (Join-Path $scriptDir "plot_kelvin_helmholtz_growth.py")
    }
    if (-not $KepsOnly -and -not $LesOnly) {
      & $python (Join-Path $scriptDir "plot_kelvin_helmholtz_snapshots.py") "pure"
    }
    if (-not $PureOnly -and -not $LesOnly) {
      & $python (Join-Path $scriptDir "plot_kelvin_helmholtz_snapshots.py") "keps"
    }
    if (-not $PureOnly -and -not $KepsOnly) {
      & $python (Join-Path $scriptDir "plot_kelvin_helmholtz_snapshots.py") "les"
    }
  }
}
finally { Pop-Location }

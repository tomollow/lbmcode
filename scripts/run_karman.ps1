param(
  [switch]$SkipPlot,
  [switch]$PureOnly,
  [switch]$KepsOnly,
  [switch]$LesOnly,
  [switch]$DesOnly
)

$ErrorActionPreference = "Stop"

$onlyCount = @($PureOnly, $KepsOnly, $LesOnly, $DesOnly | Where-Object { $_ }).Count
if ($onlyCount -gt 1) {
  throw "Specify at most one of -PureOnly / -KepsOnly / -LesOnly / -DesOnly"
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir

$variants = @()
if (-not $KepsOnly -and -not $LesOnly -and -not $DesOnly) {
  $variants += @{ Name = "karman";      Source = "src/sec4/karman.c";      Exe = "karman.exe" }
}
if (-not $PureOnly -and -not $LesOnly -and -not $DesOnly) {
  $variants += @{ Name = "karman_keps"; Source = "src/sec4/karman_keps.c"; Exe = "karman_keps.exe" }
}
if (-not $PureOnly -and -not $KepsOnly -and -not $DesOnly) {
  $variants += @{ Name = "karman_les";  Source = "src/sec4/karman_les.c";  Exe = "karman_les.exe" }
}
if (-not $PureOnly -and -not $KepsOnly -and -not $LesOnly) {
  $variants += @{ Name = "karman_des";  Source = "src/sec4/karman_des.c";  Exe = "karman_des.exe" }
}

Push-Location $repoRoot
try {
  foreach ($v in $variants) {
    & (Join-Path $scriptDir "build_one.cmd") $v.Source
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    $exePath = Join-Path $repoRoot "build/bin/$($v.Exe)"
    if (-not (Test-Path $exePath)) { throw "Executable was not generated: $exePath" }

    $runDir = Join-Path $repoRoot "outputs/sec4/$($v.Name)"
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null

    Write-Host "Running $($v.Exe) in $runDir"
    Push-Location $runDir
    try {
      & $exePath
      if ($LASTEXITCODE -ne 0) { throw "$($v.Exe) exited with code $LASTEXITCODE" }
    }
    finally { Pop-Location }
  }

  if (-not $SkipPlot) {
    $venvPython = Join-Path $repoRoot ".venv/Scripts/python.exe"
    $python = if (Test-Path $venvPython) { $venvPython } else { "python" }
    Write-Host "Regenerating plots with $python"
    if (-not $KepsOnly -and -not $LesOnly -and -not $DesOnly) {
      & $python (Join-Path $scriptDir "plot_karman_snapshots.py") "pure"
    }
    if (-not $PureOnly -and -not $LesOnly -and -not $DesOnly) {
      & $python (Join-Path $scriptDir "plot_karman_snapshots.py") "keps"
    }
    if (-not $PureOnly -and -not $KepsOnly -and -not $DesOnly) {
      & $python (Join-Path $scriptDir "plot_karman_snapshots.py") "les"
    }
    if (-not $PureOnly -and -not $KepsOnly -and -not $LesOnly) {
      & $python (Join-Path $scriptDir "plot_karman_snapshots.py") "des"
    }
    if ($onlyCount -eq 0) {
      & $python (Join-Path $scriptDir "plot_karman_spectrum.py")
    }
  }
}
finally { Pop-Location }

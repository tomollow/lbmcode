param(
  [switch]$SkipPlot,
  [switch]$KepsOnly,
  [switch]$LesOnly
)

$ErrorActionPreference = "Stop"

if ($KepsOnly -and $LesOnly) {
  throw "Cannot specify both -KepsOnly and -LesOnly"
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir

# kelbm has no separate "pure" variant — kelbm.c itself is the k-eps version,
# and kelbm_les.c is the Smagorinsky version. By default we run both.
$variants = @()
if (-not $LesOnly) {
  $variants += @{ Name = "kelbm";     Source = "src/sec4/kelbm.c";     Exe = "kelbm.exe" }
}
if (-not $KepsOnly) {
  $variants += @{ Name = "kelbm_les"; Source = "src/sec4/kelbm_les.c"; Exe = "kelbm_les.exe" }
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
    if (-not $LesOnly) {
      # k-eps plots read outputs/sec4/kelbm/kelbm_output.csv
      & $python (Join-Path $scriptDir "plot_kelbm_compare_theory.py")
      & $python (Join-Path $scriptDir "plot_kelbm_contour.py") "keps"
      & $python (Join-Path $scriptDir "plot_kelbm_centerline.py")
      & $python (Join-Path $scriptDir "plot_kelbm_output.py")
    }
    if (-not $KepsOnly) {
      # LES contour reads outputs/sec4/kelbm_les/kelbm_les_output.csv
      & $python (Join-Path $scriptDir "plot_kelbm_contour.py") "les"
    }
  }
}
finally {
  Pop-Location
}

param(
  [switch]$SkipPlot,
  [switch]$PureOnly,
  [switch]$KepsOnly
)

$ErrorActionPreference = "Stop"

if ($PureOnly -and $KepsOnly) {
  throw "Cannot specify both -PureOnly and -KepsOnly"
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir

$variants = @()
if (-not $KepsOnly) {
  $variants += @{ Name = "cavity";      Source = "src/sec4/cavity.c";      Exe = "cavity.exe" }
}
if (-not $PureOnly) {
  $variants += @{ Name = "cavity_keps"; Source = "src/sec4/cavity_keps.c"; Exe = "cavity_keps.exe" }
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
      & $python (Join-Path $scriptDir "plot_cavity_streamlines.py")
      & $python (Join-Path $scriptDir "plot_cavity_centerline.py")
    }
  }
}
finally { Pop-Location }

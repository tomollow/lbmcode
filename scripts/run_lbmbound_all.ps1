param(
  [string]$Source = "src/sec2/lbmbound.c",
  [string]$OutputRoot = "outputs/sec2/lbmbound"
)

$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir

$cases = @(
  @{ Flag = 1; Name = "equilibrium" },
  @{ Flag = 2; Name = "on_grid_bounce_back" },
  @{ Flag = 3; Name = "inamuro_no_slip" },
  @{ Flag = 4; Name = "zou_non_equilibrium" },
  @{ Flag = 5; Name = "half_way_bounce_back" },
  @{ Flag = 6; Name = "interpolated_linear" },
  @{ Flag = 7; Name = "interpolated_quadratic" }
)

Push-Location $repoRoot
try {
  & (Join-Path $scriptDir "build_one.cmd") $Source
  if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
  }

  $exePath = Join-Path $repoRoot "build/bin/lbmbound.exe"
  if (-not (Test-Path $exePath)) {
    throw "Executable was not generated: $exePath"
  }

  foreach ($case in $cases) {
    $runDir = Join-Path $repoRoot (Join-Path $OutputRoot $case.Name)
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null

    Write-Host "Running lbmbound.exe with flag $($case.Flag) in $runDir"

    Push-Location $runDir
    try {
      $case.Flag.ToString() | & $exePath *> run.log
      if ($LASTEXITCODE -ne 0) {
        exit $LASTEXITCODE
      }
    }
    finally {
      Pop-Location
    }
  }
}
finally {
  Pop-Location
}
param(
  [switch]$SkipPlot,
  [switch]$SkipRuns
)

$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot  = Split-Path -Parent $scriptDir
$source    = Join-Path $repoRoot "src/sec4/lbmcm.c"

# (Flag, Re, Subdirectory) tuples driving the comparison shown in docs/sec4/lbmcm.md.
$cases = @(
  @{ Flag = 1; Re =  100; Dir = "srt_re100"  },
  @{ Flag = 2; Re =  100; Dir = "mrt_re100"  },
  @{ Flag = 3; Re =  100; Dir = "cm_re100"   },
  @{ Flag = 2; Re = 1000; Dir = "mrt_re1000" },
  @{ Flag = 3; Re = 1000; Dir = "cm_re1000"  },
  @{ Flag = 3; Re = 5000; Dir = "cm_re5000"  }
)

function Set-LbmcmCase {
  param([string]$Content, [int]$Flag, [int]$Re)
  $result = $Content
  # Comment out every flag/re assignment line, then uncomment the requested one.
  $result = $result -replace '(?m)^  (flag = \d+;)',     '//  $1'
  $result = $result -replace '(?m)^  (re =\s*\d+;)',     '//  $1'
  $result = $result -replace "(?m)^//\s+(flag = $Flag;)", '  $1'
  $result = $result -replace "(?m)^//\s+(re =\s*$Re;)",   '  $1'
  return $result
}

Push-Location $repoRoot
try {
  if (-not $SkipRuns) {
    $original = Get-Content $source -Raw
    try {
      foreach ($c in $cases) {
        Write-Host "=== flag=$($c.Flag), re=$($c.Re), dir=$($c.Dir) ==="

        $patched = Set-LbmcmCase -Content $original -Flag $c.Flag -Re $c.Re
        Set-Content -Path $source -Value $patched -NoNewline -Encoding ascii

        & (Join-Path $scriptDir "build_one.cmd") "src/sec4/lbmcm.c"
        if ($LASTEXITCODE -ne 0) { throw "Build failed for flag=$($c.Flag), re=$($c.Re)" }

        $exe = Join-Path $repoRoot "build/bin/lbmcm.exe"
        if (-not (Test-Path $exe)) { throw "Executable not found: $exe" }

        $runDir = Join-Path $repoRoot "outputs/sec4/lbmcm/$($c.Dir)"
        New-Item -ItemType Directory -Path $runDir -Force | Out-Null

        Push-Location $runDir
        try {
          & $exe | Out-File -FilePath run.log -Encoding utf8
          if ($LASTEXITCODE -ne 0) { throw "lbmcm.exe exited $LASTEXITCODE in $runDir" }
        }
        finally { Pop-Location }
      }
    }
    finally {
      # Always restore the original source so working tree is clean.
      Set-Content -Path $source -Value $original -NoNewline -Encoding ascii
      Write-Host "Restored $source"
    }
  }

  if (-not $SkipPlot) {
    $venvPython = Join-Path $repoRoot ".venv/Scripts/python.exe"
    $python = if (Test-Path $venvPython) { $venvPython } else { "python" }
    & $python (Join-Path $scriptDir "plot_lbmcm_distribution.py")
    & $python (Join-Path $scriptDir "plot_lbmcm_compare.py")
  }
}
finally { Pop-Location }

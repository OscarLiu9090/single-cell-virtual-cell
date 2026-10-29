param([Parameter(Mandatory=$true)][string]$Config, [switch]$Preflight)
$ErrorActionPreference = 'Stop'
$skillRoot = Split-Path -Parent $PSScriptRoot
$runtimeFile = Join-Path $skillRoot 'references\environment.local.json'
if (!(Test-Path -LiteralPath $runtimeFile)) { throw 'Configure references/environment.local.json with the Rscript executable path.' }
$runtime = Get-Content -LiteralPath $runtimeFile -Raw -Encoding UTF8 | ConvertFrom-Json
$entry = Join-Path $PSScriptRoot 'run_perturbation.R'
$resolvedConfig = (Resolve-Path -LiteralPath $Config).Path
if ($Preflight) { & $runtime.rscript --vanilla $entry --config $resolvedConfig --preflight }
else { & $runtime.rscript --vanilla $entry --config $resolvedConfig }
exit $LASTEXITCODE

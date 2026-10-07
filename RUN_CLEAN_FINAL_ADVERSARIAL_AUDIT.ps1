param(
    [string]$Rscript = "",
    [switch]$NoParallel,
    [switch]$SkipReport
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Invoke-NativeChecked {
    param(
        [Parameter(Mandatory=$true)][string]$Exe,
        [Parameter(ValueFromRemainingArguments=$true)][string[]]$Args
    )
    & $Exe @Args
    if ($LASTEXITCODE -ne 0) {
        throw "Native command failed ($LASTEXITCODE): $Exe $($Args -join ' ')"
    }
}

$root = (& git rev-parse --show-toplevel).Trim()
if ($LASTEXITCODE -ne 0 -or -not $root) { throw "Not inside the NIST git repository." }
Set-Location $root

$dirty = @(git status --porcelain)
if ($dirty.Count -ne 0) {
    Write-Host "Working tree must be clean before the canonical release run:" -ForegroundColor Red
    $dirty | ForEach-Object { Write-Host $_ }
    throw "Clean-start gate failed."
}

if (-not $Rscript) {
    if ($env:RSCRIPT -and (Test-Path $env:RSCRIPT)) {
        $Rscript = $env:RSCRIPT
    } else {
        $rDir = Get-ChildItem "C:\Program Files\R" -Directory -ErrorAction SilentlyContinue |
            Sort-Object Name -Descending |
            Select-Object -First 1
        if (-not $rDir) { throw "R installation not found. Pass -Rscript <path>." }
        $Rscript = Join-Path $rDir.FullName "bin\Rscript.exe"
    }
}
if (-not (Test-Path $Rscript)) { throw "Rscript not found: $Rscript" }

$pinnedCommit = (& git rev-parse HEAD).Trim()
if ($pinnedCommit -notmatch '^[0-9a-fA-F]{40}$') { throw "Could not resolve a full pinned commit SHA." }

Write-Host "Pinned clean input commit: $pinnedCommit"
Write-Host "Verifying frozen NIST source SHA-256 provenance..."
Invoke-NativeChecked $Rscript "R/verify_source_provenance.R" "--check"

# Preserve the clean-start attestation for the immutable-style snapshot step.
$env:NIST_AUDIT_INPUT_COMMIT = $pinnedCommit
$env:NIST_AUDIT_CLEAN_START = "TRUE"

$parallel = if ($NoParallel) { "false" } else { "true" }
$render = if ($SkipReport) { "false" } else { "true" }

Write-Host "Running canonical strict final adversarial audit..."
Invoke-NativeChecked $Rscript `
    "R/render_statistical_audit.R" `
    "--bootstrap-n" "5000" `
    "--null-n" "5000" `
    "--calibration-n" "10000" `
    "--injection-n" "2000" `
    "--seed" "20260517" `
    "--parallel" $parallel `
    "--force" `
    "--strict" `
    "--render-report" $render

$jsonPath = Join-Path $root "final_adversarial_audit.json"
if (-not (Test-Path $jsonPath)) { throw "Canonical run did not produce final_adversarial_audit.json." }
$result = Get-Content $jsonPath -Raw | ConvertFrom-Json

if ($result.git_commit -ne $pinnedCommit) {
    throw "Machine result commit mismatch: expected $pinnedCommit, got $($result.git_commit)"
}
if ([int]$result.canonical_seed -ne 20260517) { throw "Canonical seed mismatch." }
if ([int]$result.budgets.bootstrap_n -lt 5000 -or
    [int]$result.budgets.null_n -lt 5000 -or
    [int]$result.budgets.calibration_n -lt 10000 -or
    [int]$result.budgets.injection_n -lt 2000) {
    throw "Machine result does not satisfy the frozen Monte-Carlo budgets."
}

# A canonical run necessarily writes declared output artifacts. It must not
# alter source code, input data, protocol configuration, or other undeclared
# files. Reject any such mutation before snapshotting.
$changed = @()
$changed += @(git diff --name-only)
$changed += @(git ls-files --others --exclude-standard)
$changed = @($changed | Where-Object { $_ } | Sort-Object -Unique)

$allowedPatterns = @(
    '^FINAL_ADVERSARIAL_AUDIT\.md$',
    '^REPRODUCE_FINAL_ADVERSARIAL_AUDIT\.txt$',
    '^final_adversarial_audit\.json$',
    '^tables_r/statistical_audit/',
    '^figures_r/statistical_audit/',
    '^outputs_r/statistical_audit/',
    '^reports/rendered/',
    '^results/final_adversarial_audit_v1/?'
)

$unexpected = @()
foreach ($p in $changed) {
    $norm = $p -replace '\\','/'
    $ok = $false
    foreach ($pattern in $allowedPatterns) {
        if ($norm -match $pattern) { $ok = $true; break }
    }
    if (-not $ok) { $unexpected += $p }
}
if ($unexpected.Count -gt 0) {
    Write-Host "Unexpected files changed during canonical run:" -ForegroundColor Red
    $unexpected | ForEach-Object { Write-Host $_ }
    throw "Release source-integrity gate failed."
}

Write-Host "Creating verdict-preserving immutable-style release snapshot..."
Invoke-NativeChecked $Rscript "R/freeze_release_snapshot.R" "--input-commit" $pinnedCommit

$classification = [string]$result.classifications.Fe
$tagSuffix = if ($classification -like 'PASS*') { 'pass' } elseif ($classification -eq 'CONDITIONAL') { 'conditional' } else { 'fail' }
$tag = "nist-final-adversarial-audit-v1-$tagSuffix"

Write-Host ""
Write-Host "Canonical run completed from clean pinned input commit." -ForegroundColor Green
Write-Host "Fe II classification: $classification"
Write-Host "Co II classification: $($result.classifications.Co)"
Write-Host ""
Write-Host "Review the generated diff, then commit the release outputs. After that commit, create the immutable tag:" -ForegroundColor Yellow
Write-Host "  git add FINAL_ADVERSARIAL_AUDIT.md REPRODUCE_FINAL_ADVERSARIAL_AUDIT.txt final_adversarial_audit.json tables_r/statistical_audit figures_r/statistical_audit outputs_r/statistical_audit reports/rendered results"
Write-Host "  git commit -m `"Freeze final NIST adversarial audit ($tagSuffix)`""
Write-Host "  git tag -a $tag -m `"Final NIST adversarial audit: $classification`""
Write-Host "  git push origin HEAD"
Write-Host "  git push origin $tag"
Write-Host ""
Write-Host "Current status:"
git status --short

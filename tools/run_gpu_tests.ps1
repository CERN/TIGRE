<#
.SYNOPSIS
Run TIGRE's GPU and MATLAB tests on this machine, for a branch, a pull request
or a commit.

.DESCRIPTION
GitHub's runners have no GPU, and a self-hosted runner on a public repository
could run code from fork pull requests. So the GPU tests are run by hand with
this script instead.

It checks the code out into a separate folder (a git worktree), so your own
working tree is never touched, then:
  1. builds and installs TIGRE in a fresh virtual environment, only for this
     machine's GPUs (TIGRE_CUDA_ARCH=native, a much faster build);
  2. runs the Python tests (GPU tests included; -Extended adds the stress tests);
  3. compiles the MATLAB MEX files and runs every MATLAB/Tests/test_*.m
     (skip with -NoMatlab).
At the end it prints a short summary, also saved as summary.md next to the logs,
ready to paste into the pull request.

Read what a pull request changes before testing it: the script builds and runs
that code on this machine.

.EXAMPLE
From cmd or the Anaconda prompt, inside the repository:
  powershell -ExecutionPolicy Bypass -File tools\run_gpu_tests.ps1 -Ref 944

.EXAMPLE
.\tools\run_gpu_tests.ps1                    # master
.\tools\run_gpu_tests.ps1 -Ref 944           # pull request #944
.\tools\run_gpu_tests.ps1 -Ref fix/my-branch -Extended
.\tools\run_gpu_tests.ps1 -Ref 944 -NoMatlab
#>
param(
    # Branch name, pull request number or commit. Default: master.
    [string]$Ref = "master",
    # Also run the tests marked `extended` (stress tests, much slower).
    [switch]$Extended,
    # Skip the MATLAB compile and tests.
    [switch]$NoMatlab,
    # Python to build the virtual environment with.
    [string]$Python = "python",
    # Where the checkout, environment and logs go.
    [string]$WorkDir = (Join-Path $env:TEMP "tigre-gpu-tests"),
    # Keep the checkout afterwards (it is removed by default; logs are kept).
    [switch]$Keep,
    # The TIGRE repository to fetch from. Default: the one the current folder
    # is in, else the one this script is in.
    [string]$Repo = ""
)

$ErrorActionPreference = "Continue"
if (-not $Repo) {
    $Repo = git rev-parse --show-toplevel 2>$null
    if (-not $Repo) { $Repo = git -C $PSScriptRoot rev-parse --show-toplevel 2>$null }
}
if (-not $Repo) {
    Write-Host "ERROR: run this from inside the TIGRE repository, or pass -Repo D:\path\to\TIGRE" -ForegroundColor Red
    exit 1
}
$repo = "$Repo".Trim()
$src = Join-Path $WorkDir "src"
$stamp = Get-Date -Format "yyyy-MM-dd_HHmm"
$logs = Join-Path $WorkDir "logs\$stamp"
New-Item -ItemType Directory -Force -Path $logs | Out-Null

function Stop-WithError($message) {
    Write-Host "ERROR: $message" -ForegroundColor Red
    exit 1
}

# 1. Fetch the code to test.
if ($Ref -match '^\d+$') {
    $fetchRef = "pull/$Ref/head"
    $label = "PR #$Ref"
} else {
    $fetchRef = $Ref
    $label = $Ref
}
Write-Host "== Fetching $label"
git -C $repo fetch --quiet origin $fetchRef
if ($LASTEXITCODE -ne 0) { Stop-WithError "could not fetch '$fetchRef' from origin" }
if (Test-Path $src) {
    git -C $repo worktree remove --force $src 2>$null
    if (Test-Path $src) { Remove-Item -Recurse -Force $src }
    git -C $repo worktree prune
}
git -C $repo worktree add --quiet --detach $src FETCH_HEAD
if ($LASTEXITCODE -ne 0) { Stop-WithError "could not create the worktree in $src" }
$commit = (git -C $src rev-parse --short HEAD).Trim()
Write-Host "   commit $commit, checked out in $src"

# 2. Build and install TIGRE.
Write-Host "== Building TIGRE (TIGRE_CUDA_ARCH=native)"
$env:TIGRE_CUDA_ARCH = "native"
$env:MPLBACKEND = "Agg"
$venv = Join-Path $WorkDir "venv"
if (Test-Path $venv) { Remove-Item -Recurse -Force $venv }
& $Python -m venv $venv
if ($LASTEXITCODE -ne 0) { Stop-WithError "could not create a virtual environment with '$Python'" }
$py = Join-Path $venv "Scripts\python.exe"
& $py -m pip install --quiet --upgrade pip *> (Join-Path $logs "pip.log")
& $py -m pip install $src pytest *>> (Join-Path $logs "pip.log")
if ($LASTEXITCODE -ne 0) { Stop-WithError "the build failed; see $logs\pip.log" }
$versions = (& $py -c "import tigre; tigre.show_versions()" 2>&1 | Out-String).Trim()

# 3. Python tests. Run from the checkout so pytest uses its pyproject.toml settings.
Write-Host "== Python tests"
Push-Location $src
if ($Extended) {
    & $py -m pytest -ra -m "extended or not extended" Python/tests 2>&1 | Tee-Object (Join-Path $logs "pytest.log")
} else {
    & $py -m pytest -ra Python/tests 2>&1 | Tee-Object (Join-Path $logs "pytest.log")
}
$pyCode = $LASTEXITCODE
Pop-Location
$pyLast = (Get-Content (Join-Path $logs "pytest.log") | Where-Object { $_ -match '(passed|failed|error)' } | Select-Object -Last 1)
$pyResult = if ($pyCode -eq 0) { "PASS" } else { "FAIL" }

# 4. MATLAB: compile, then each test function. Every test runs; failures are counted.
$mlResult = "skipped"
if (-not $NoMatlab) {
    Write-Host "== MATLAB compile and tests"
    $matlabScript = @"
cd('$($src -replace "'", "''")\MATLAB');
setenv('TIGRE_CUDA_ARCH','native');
try, Compile; catch e, disp(getReport(e)); exit(2); end
InitTIGRE; addpath('Tests');
t = dir('Tests/test_*.m'); failed = {};
for k = 1:numel(t)
    name = t(k).name(1:end-2);
    try, feval(name); catch e, failed{end+1} = name; fprintf('FAIL %s\n%s\n', name, getReport(e)); end
end
fprintf('MATLAB tests: %d run, %d failed\n', numel(t), numel(failed));
exit(double(~isempty(failed)));
"@
    $mFile = Join-Path $WorkDir "run_matlab_tests.m"
    Set-Content -Path $mFile -Value $matlabScript -Encoding ASCII
    matlab -batch "run('$($mFile -replace "'", "''")')" 2>&1 | Tee-Object (Join-Path $logs "matlab.log")
    switch ($LASTEXITCODE) {
        0 { $mlResult = "PASS" }
        2 { $mlResult = "FAIL (compile)" }
        default { $mlResult = "FAIL" }
    }
    $mlLast = (Get-Content (Join-Path $logs "matlab.log") | Where-Object { $_ -match '^MATLAB tests:' } | Select-Object -Last 1)
    if ($mlLast) { $mlResult = "$mlResult ($mlLast)" }
}

# 5. Summary.
$mode = if ($Extended) { "core + extended" } else { "core" }
$summary = @"
**GPU tests** for $label (commit $commit), $(Get-Date -Format "yyyy-MM-dd HH:mm"), run with ``tools/run_gpu_tests.ps1``:

- Python ($mode): **$pyResult**, $pyLast
- MATLAB: **$mlResult**

<details><summary>Machine</summary>

``````
$versions
``````
</details>
"@
Set-Content -Path (Join-Path $logs "summary.md") -Value $summary -Encoding UTF8
Write-Host ""
Write-Host $summary
Write-Host ""
Write-Host "Logs and summary.md: $logs"

if (-not $Keep) {
    git -C $repo worktree remove --force $src 2>$null
    git -C $repo worktree prune
}
if ($pyResult -ne "PASS" -or ($mlResult -ne "skipped" -and -not $mlResult.StartsWith("PASS"))) { exit 1 }
exit 0

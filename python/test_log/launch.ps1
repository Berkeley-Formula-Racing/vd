[CmdletBinding()]
param(
    [string]$DataDir,
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$ForwardedArguments
)

$callerDirectory = (Get-Location).ProviderPath
$resolvedDataDir = $null
if ($PSBoundParameters.ContainsKey("DataDir")) {
    if ([System.IO.Path]::IsPathRooted($DataDir)) {
        $resolvedDataDir = [System.IO.Path]::GetFullPath($DataDir)
    }
    else {
        $resolvedDataDir = [System.IO.Path]::GetFullPath(
            (Join-Path -Path $callerDirectory -ChildPath $DataDir)
        )
    }
}

$candidates = @()
$pyLauncher = Get-Command -Name "py" -CommandType Application -ErrorAction SilentlyContinue |
    Select-Object -First 1
if ($null -ne $pyLauncher) {
    $candidates += [pscustomobject]@{
        Command = $pyLauncher.Source
        PrefixArguments = @("-3.12")
        Label = "py -3.12"
    }
}

$pythonCommand = Get-Command -Name "python" -CommandType Application -ErrorAction SilentlyContinue |
    Select-Object -First 1
if ($null -ne $pythonCommand) {
    $candidates += [pscustomobject]@{
        Command = $pythonCommand.Source
        PrefixArguments = @()
        Label = "python"
    }
}

if (-not [string]::IsNullOrWhiteSpace($env:USERPROFILE)) {
    $bundledRuntime = Join-Path -Path $env:USERPROFILE -ChildPath ".cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
    if (Test-Path -LiteralPath $bundledRuntime -PathType Leaf) {
        $candidates += [pscustomobject]@{
            Command = $bundledRuntime
            PrefixArguments = @()
            Label = "bundled Codex Python"
        }
    }
}

$preflight = @"
import sqlite3
import sys

if sys.version_info < (3, 12):
    sys.stderr.write('Python 3.12 or newer is required; found ' + sys.version.split()[0] + '\n')
    raise SystemExit(10)

import tkinter as tk

root = None
try:
    root = tk.Tk()
    root.withdraw()
    root.update_idletasks()
finally:
    if root is not None:
        try:
            root.destroy()
        except tk.TclError:
            pass

print('Python, sqlite3, and a withdrawn Tk window are available.')
"@

$diagnostics = New-Object 'System.Collections.Generic.List[string]'
$selected = $null
foreach ($candidate in $candidates) {
    $probeArguments = @($candidate.PrefixArguments) + @("-c", $preflight)
    $probeOutput = & $candidate.Command @probeArguments 2>&1
    $probeExitCode = $LASTEXITCODE
    if ($probeExitCode -eq 0) {
        $selected = $candidate
        break
    }

    $detail = (($probeOutput | ForEach-Object { "$_" }) -join [Environment]::NewLine).Trim()
    if ([string]::IsNullOrWhiteSpace($detail)) {
        $detail = "The runtime exited with code $probeExitCode."
    }
    $null = $diagnostics.Add("$($candidate.Label): $detail")
}

if ($null -eq $selected) {
    [Console]::Error.WriteLine(
        "No compatible Python 3.12+ runtime with sqlite3 and working Tk support was found. Install Python 3.12 or newer with Tcl/Tk enabled, or repair the bundled runtime."
    )
    foreach ($diagnostic in $diagnostics) {
        [Console]::Error.WriteLine("  $diagnostic")
    }
    exit 1
}

$pythonArguments = @()
if ($null -ne $resolvedDataDir) {
    $pythonArguments += @("--data-dir", $resolvedDataDir)
}
if ($null -ne $ForwardedArguments) {
    $pythonArguments += $ForwardedArguments
}

$originalDirectory = (Get-Location).ProviderPath
$launchExitCode = 1
try {
    Set-Location -LiteralPath $PSScriptRoot -ErrorAction Stop
    $moduleArguments = @($selected.PrefixArguments) + @("-m", "vd_test_log") + $pythonArguments
    & $selected.Command @moduleArguments
    $launchExitCode = $LASTEXITCODE
    if ($null -eq $launchExitCode) {
        $launchExitCode = 1
    }
}
catch {
    [Console]::Error.WriteLine("Could not launch the test log with $($selected.Label): $_")
    $launchExitCode = 1
}
finally {
    if (-not [string]::IsNullOrWhiteSpace($originalDirectory)) {
        Set-Location -LiteralPath $originalDirectory -ErrorAction SilentlyContinue
    }
}

exit [int]$launchExitCode
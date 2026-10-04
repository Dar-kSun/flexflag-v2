# Docking validation (docs/plan-v0.6-docking.md), unattended. About 10 hours on 18 cores.
#
#   powershell -ExecutionPolicy Bypass -File scripts\overnight_docking.ps1
#
# Keeps the PC awake while running. Resumable: proteins already in
# results\docking\runs.csv are skipped, so re-running after an interruption continues
# where it stopped. Progress: data\cache\docking.log. When done it writes
# results\docking\summary.md.

$ErrorActionPreference = "Continue"
$Repo = Split-Path -Parent $PSScriptRoot
Set-Location $Repo
$env:PYTHONUTF8 = "1"
$Py = Join-Path $Repo ".venv\Scripts\python.exe"
$Log = Join-Path $Repo "data\cache\docking.log"

Add-Type -Namespace Win32 -Name Power -MemberDefinition @'
[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);
'@
[Win32.Power]::SetThreadExecutionState([uint32]"0x80000001") | Out-Null  # stay awake

"docking run started $(Get-Date)" | Out-File $Log -Append -Encoding utf8
# Progress lines only; RDKit's per-molecule warnings are dropped from the log.
& $Py -u scripts\14_docking.py --workers 18 2>&1 |
    ForEach-Object { "$_" } |
    Where-Object { $_ -notmatch "tagged as 2D|UFFTYPER" } |
    Tee-Object -Variable lines | Out-Host
$lines | Out-File $Log -Append -Encoding utf8
& $Py -u scripts\14_docking.py --summarise 2>&1 | ForEach-Object { "$_" } |
    Out-File $Log -Append -Encoding utf8
"docking run finished $(Get-Date)" | Out-File $Log -Append -Encoding utf8

[Win32.Power]::SetThreadExecutionState([uint32]"0x80000000") | Out-Null  # back to normal

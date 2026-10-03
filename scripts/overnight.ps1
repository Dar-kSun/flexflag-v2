# Overnight run: the analyses pre-declared in docs/plan-v0.5-overnight.md.
#
#   powershell -ExecutionPolicy Bypass -File scripts\overnight.ps1
#
# Runs unattended and keeps the PC awake while it runs. Each step is independent:
# a failure is logged and the next step still runs. Everything goes to
# data\cache\overnight.log; a per-step status table goes to
# data\cache\overnight_status.txt. Safe to re-run: every step resumes from its cache.

$ErrorActionPreference = "Continue"
$Repo = Split-Path -Parent $PSScriptRoot
Set-Location $Repo
$env:PYTHONUTF8 = "1"
$Py = Join-Path $Repo ".venv\Scripts\python.exe"
$Log = Join-Path $Repo "data\cache\overnight.log"
$Status = Join-Path $Repo "data\cache\overnight_status.txt"
New-Item -ItemType Directory -Force (Join-Path $Repo "data\cache") | Out-Null

# Keep the system (not the display) awake until this script exits.
Add-Type -Namespace Win32 -Name Power -MemberDefinition @'
[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);
'@
[Win32.Power]::SetThreadExecutionState([uint32]"0x80000001") | Out-Null  # CONTINUOUS | SYSTEM_REQUIRED

function Step($Name, [scriptblock]$Body) {
    $t0 = Get-Date
    "`n===== $Name  (started $t0) =====" | Tee-Object -FilePath $Log -Append
    & $Body 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $Log -Append
    $code = $LASTEXITCODE
    $mins = [math]::Round(((Get-Date) - $t0).TotalMinutes, 1)
    $line = "{0,-34} exit={1,-4} {2} min" -f $Name, $code, $mins
    $line | Tee-Object -FilePath $Status -Append
}

"overnight run started $(Get-Date)" | Out-File $Status -Encoding utf8

Step "1 install torch (CUDA) + fair-esm" {
    & $Py -c "import torch, esm" 2>$null
    if ($LASTEXITCODE -ne 0) {
        & $Py -m pip install --progress-bar off torch --index-url https://download.pytorch.org/whl/cu126
        & $Py -m pip install --progress-bar off fair-esm
    }
    & $Py -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"
}
Step "2 ESM-2 embeddings (GPU)"        { & $Py -u scripts\11_esm_embed.py }
Step "3 ESM-2 vs pocket pLDDT"         { & $Py -u scripts\12_esm_eval.py }
Step "4 many ligands per protein"      { & $Py -u scripts\13_multipair.py --workers 24 }
Step "5 tests still pass"              { & $Py -m pytest -q }

"overnight run finished $(Get-Date)" | Tee-Object -FilePath $Status -Append
[Win32.Power]::SetThreadExecutionState([uint32]"0x80000000") | Out-Null      # back to normal

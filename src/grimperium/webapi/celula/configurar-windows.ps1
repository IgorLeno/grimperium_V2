# Lado Windows da célula do Grimperium. O instalador (instalar.sh, no WSL)
# preenche a distribuição e o usuário e grava este arquivo na sua pasta do
# Windows. Rode no PowerShell do seu usuário (não precisa ser administrador):
#
#   powershell -ExecutionPolicy Bypass -File "$env:USERPROFILE\grimperium-celula-windows.ps1"
#
# 1. .wslconfig: o WSL usa todos os núcleos e 75% da memória (cópia do
#    arquivo anterior em .wslconfig.antes-grimperium, feita só uma vez).
# 2. Tarefa "Grimperium Celula" no Agendador: inicia a célula no logon deste
#    usuário, sem janela, sem limite de tempo, também na bateria.
# 3. Energia: o computador não suspende nem hiberna (tela pode desligar).
# 4. Atalho "Grimperium Célula" na Área de Trabalho: abre o painel ao vivo
#    (celula painel). Fechar a janela não para o cálculo.
#
# Depois, reinicie o Windows. Para a célula voltar sozinha depois de uma
# queda de energia sem ninguém fazer logon, ative o logon automático
# (netplwiz) e, na BIOS, "ligar ao voltar a energia" (veja o guia).

$ErrorActionPreference = 'Stop'

$Distro = '@DISTRO@'
$LinuxUser = '@LINUX_USER@'
$Iniciar = '@INICIAR@'
$Painel = '@PAINEL@'
$TaskName = 'Grimperium Celula'

if ($Distro.StartsWith('@')) {
    throw 'Este é o modelo do pacote. Use o arquivo que o instalar.sh gravou na sua pasta do Windows.'
}

# --- 1. .wslconfig ------------------------------------------------------------
$computer = Get-CimInstance Win32_ComputerSystem
$cores = [int]$computer.NumberOfLogicalProcessors
$ramGB = [math]::Floor($computer.TotalPhysicalMemory / 1GB)
$memGB = [math]::Max(2, [math]::Floor($ramGB * 0.75))
$wslconfig = Join-Path $env:USERPROFILE '.wslconfig'
$backup = "$wslconfig.antes-grimperium"
# Só a primeira execução copia: numa segunda, o .wslconfig já é o nosso e
# sobrescreveria a cópia do original.
if ((Test-Path $wslconfig) -and -not (Test-Path $backup)) {
    Copy-Item $wslconfig $backup
    Write-Host "Cópia do .wslconfig anterior: $backup"
}
$lines = @(
    '# Escrito pelo instalador da célula do Grimperium.',
    '[wsl2]',
    "processors=$cores",
    "memory=${memGB}GB"
)
Set-Content -Path $wslconfig -Value $lines -Encoding ascii
Write-Host "1/4 .wslconfig: $cores núcleos, ${memGB} GB de memória para o WSL."

# --- 2. tarefa no Agendador ---------------------------------------------------
$wsl = Join-Path $env:WINDIR 'System32\wsl.exe'
$user = "$env:USERDOMAIN\$env:USERNAME"
# PowerShell oculto só para não abrir uma janela que alguém feche sem querer.
$command = "& '$wsl' -d '$Distro' -u '$LinuxUser' -- '$Iniciar'"
$action = New-ScheduledTaskAction -Execute 'powershell.exe' `
    -Argument "-NoProfile -WindowStyle Hidden -Command `"$command`""
$trigger = New-ScheduledTaskTrigger -AtLogOn -User $user
$settings = New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries `
    -DontStopIfGoingOnBatteries -StartWhenAvailable `
    -ExecutionTimeLimit ([TimeSpan]::Zero) -MultipleInstances IgnoreNew `
    -RestartCount 3 -RestartInterval (New-TimeSpan -Minutes 5)
$principal = New-ScheduledTaskPrincipal -UserId $user -LogonType Interactive -RunLevel Limited
Register-ScheduledTask -TaskName $TaskName -Action $action -Trigger $trigger `
    -Settings $settings -Principal $principal -Force | Out-Null
Write-Host "2/4 Tarefa '$TaskName' criada: inicia a célula no logon de $user."

# --- 3. energia ---------------------------------------------------------------
foreach ($setting in 'standby-timeout-ac', 'standby-timeout-dc', 'hibernate-timeout-ac', 'hibernate-timeout-dc') {
    powercfg /change $setting 0
    if ($LASTEXITCODE -ne 0) { Write-Warning "powercfg /change $setting 0 falhou" }
}
Write-Host '3/4 Energia: sem suspensão e sem hibernação.'

# --- 4. atalho do painel ------------------------------------------------------
# GetFolderPath acha a Área de Trabalho também quando o OneDrive a move.
$desktop = [Environment]::GetFolderPath('Desktop')
$shortcutPath = Join-Path $desktop 'Grimperium Célula.lnk'
$shortcut = (New-Object -ComObject WScript.Shell).CreateShortcut($shortcutPath)
$shortcut.TargetPath = $wsl
$shortcut.Arguments = "-d `"$Distro`" -u `"$LinuxUser`" -- `"$Painel`""
$shortcut.Description = 'Painel ao vivo da célula do Grimperium (fechar não para o cálculo)'
$shortcut.IconLocation = "$wsl,0"
$shortcut.Save()
Write-Host "4/4 Atalho do painel na Área de Trabalho: $shortcutPath"

Write-Host ''
Write-Host 'Pronto. Reinicie o Windows: a célula inicia sozinha depois do logon.'
Write-Host "Para iniciar agora sem reiniciar: wsl --shutdown; Start-ScheduledTask -TaskName '$TaskName'"

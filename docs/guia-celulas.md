# Guia: instalar uma célula em outro computador

Uma **célula** é um PC com Windows + WSL2 (Ubuntu) que calcula moléculas da
fila de treinamento para esta máquina, a **central**. A célula pede trabalho
pela rede privada do Tailscale; a central nunca entra na célula. Manutenção e
reinício são feitos pelo TeamViewer.

Você precisa de: acesso de administrador ao Windows da célula (só para
instalar o WSL e o Tailscale), internet nos dois computadores e uns 2 GB
livres no Ubuntu da célula.

## 1. Tailscale nos dois computadores

1. Instale o Tailscale (<https://tailscale.com/download>) na central e no
   Windows da célula e entre com **a mesma conta** nos dois.
2. Na central, anote o IP do Tailscale:

   ```bash
   tailscale ip -4
   ```

   Algo como `100.101.102.103`. O endereço da central para as células é
   `http://100.101.102.103:8100`.
3. (Recomendado) No painel do Tailscale, limite a célula à porta 8100 da
   central (ACL).

## 2. Subir a API das células na central

Ela é separada do dashboard e escuta só no IP do Tailscale. Num terminal da
central, na pasta do projeto:

```bash
poetry run uvicorn grimperium.webapi.cells:create_cell_app --factory --host 100.101.102.103 --port 8100
```

(troque pelo seu IP). Deixe esse terminal aberto. Se trocar CREST, xtb ou
MOPAC da central, reinicie esta API.

## 3. Criar a célula e baixar o pacote

1. Abra o dashboard: **Calcular PM7 → Células**.
2. Confira o **endereço da central** (o dashboard sugere o IP do Tailscale
   com a porta 8100), dê um nome ao computador e clique em **Criar célula e
   baixar pacote**. O navegador baixa `celula-<nome>-<id>.zip`.

**O zip contém o token da célula**: quem tiver o arquivo pode pedir moléculas
em nome dela. Leve-o para a célula por pendrive ou pelo TeamViewer, nunca por
e-mail ou pasta compartilhada, e apague-o depois de instalar.

Como o token: a central guarda só o hash dele, então o token existe apenas
dentro do zip. Por isso **cada pacote gera um token novo** e invalida o
anterior. Perdeu o zip ou precisa reinstalar? Clique em **Baixar pacote
novo** na linha da célula e reinstale com ele (a instalação antiga para com
"token recusado" até receber o pacote novo). Vazou? **Revogue** a célula e
crie outra.

## 4. Instalar na célula

No Windows da célula:

1. Instale o WSL com Ubuntu (PowerShell como administrador, uma vez):
   `wsl --install -d Ubuntu`, reinicie e crie o usuário do Ubuntu.
2. Clique com o botão direito no zip → **Extrair tudo** (em Downloads).
3. Abra o **Ubuntu** e rode (troque `SEU_USUARIO` e o nome da pasta):

   ```bash
   bash /mnt/c/Users/SEU_USUARIO/Downloads/celula-nome-id/instalar.sh
   ```

   Ele instala, sem sudo: o micromamba, o ambiente `grimperium-quimica`
   (CREST, xtb e MOPAC iguais aos da central, pelo `conda-lock.yml`), um
   Python separado para o programa `celula` e a pasta `~/grimperium-celula/`
   (nunca em `/mnt/c`, com `config.json` em `chmod 600`). No fim ele confere
   se os programas batem com os da central e se a central responde.
4. O instalador mostra um comando do PowerShell. Rode-o **no PowerShell do
   seu usuário** (não precisa ser administrador). Ele:
   - escreve `%USERPROFILE%\.wslconfig` (todos os núcleos, 75% da memória;
     o anterior fica em `.wslconfig.antes-grimperium`, copiado só na
     primeira vez, então rodar de novo não apaga o original);
   - cria a tarefa **Grimperium Celula** no Agendador, que inicia a célula no
     logon, sem janela;
   - tira a suspensão e a hibernação do plano de energia.
5. Reinicie o Windows e faça logon. A célula sobe sozinha.

Para voltar sozinha depois de uma queda de energia sem ninguém por perto:
ative o logon automático do Windows (`netplwiz`, desmarque "Os usuários
devem digitar um nome…") e, na BIOS, a opção de ligar quando a energia
volta ("Restore on AC power loss" / "Power On").

## 5. Teste de admissão

Uma célula nova fica **em admissão** e, antes de receber a fila, calcula
sozinha 5 moléculas de referência que a central já calculou: `cbs_09698` CO,
`cbs_00187` [H][C]=C, `cbs_24671` C[N+](=O)[O-], `cbs_08982` OCCO e
`cbs_01400` CCC(=O)O. A central compara cada HOF com o valor que ela aceitou:

- as 5 dentro de **0,1 kcal/mol**: a célula passa para **Admitida** e começa
  a receber a fila (em até 1 minuto);
- uma fora da tolerância, ou uma que falhe (timeout do CREST, erro do MOPAC):
  a célula fica **Bloqueada** e o dashboard mostra o motivo, por exemplo
  `cbs_09698 (CO): HOF -48.70000 contra -48.93880 da central (Δ 0.23880 > 0.1
  kcal/mol)`.

Os resultados das referências não entram na fila nem nas geometrias da
central. No dashboard (Calcular PM7 → Células), enquanto o teste roda, a
linha mostra quantas referências já voltaram. Se a célula reiniciar no meio,
ela recalcula só as que faltavam.

Depois de um bloqueio, corrija a causa (em geral, reinstale com um pacote
novo) e clique em **Refazer admissão**: os resultados anteriores são
descartados e as 5 são calculadas de novo. **Admitir** e **Bloquear** seguem
disponíveis para decidir à mão.

Se o dashboard mostrar **Versão divergente**, os programas da célula não são
os da central: rode o instalador de novo com um pacote novo.

## 6. Acompanhar na célula

No Ubuntu da célula:

```bash
celula status
```

Mostra se está rodando, a fase, a admissão, os workers ocupados, o último
contato com a central e os resultados à espera de envio. Log completo em
`~/grimperium-celula/celula.log`. O `status` funciona mesmo sem rede.

## 7. Reiniciar pelo TeamViewer

- **A célula travou ou parou**: no PowerShell da célula,
  `wsl --shutdown` e depois `Start-ScheduledTask -TaskName 'Grimperium Celula'`
  (ou reinicie o Windows). Ao voltar, ela se registra de novo, recomeça o
  cálculo interrompido (perda máxima de ~1 h de CREST) e envia o que já
  estava pronto.
- **Queda de rede**: nada a fazer; ela tenta de novo sozinha.
- **Computador desligado por muito tempo**: as moléculas dela continuam
  reservadas. No dashboard, **Devolver à fila** as libera para as outras.
- **Atualizar o programa**: baixe um pacote novo e rode o instalador de novo
  (pare a célula antes: `wsl --shutdown`).

Teste manual sem o Agendador (no Ubuntu): `~/grimperium-celula/iniciar.sh`
(Ctrl+C não alcança a célula iniciada assim; use `celula start` em primeiro
plano se quiser pará-la com Ctrl+C).

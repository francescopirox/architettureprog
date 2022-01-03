;---------------------------------------------------------
; Regressione con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf32 fss32.nasm 
;
%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati

section .bss			; Sezione contenente dati non inizializzati
	alignb 16
	stepind		resd		1


section .text			; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global prova

input		equ		8

msg	db	'stepind:',32,0
nl	db	10,0



prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; legge i parametri dal Record di Attivazione corrente
		; ------------------------------------------------------------

		; elaborazione
		
		; esempio: stampa input->stepind
		mov EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
		; [EAX]	input->x
		; [EAX + 4] input->xh
		; [EAX + 8] input->c
		; [EAX + 12] input->r
		; [EAX + 16] input->nx
		; [EAX + 20] input->d
		; [EAX + 24] input->iter
		; [EAX + 28] input->stepind
		; [EAX + 32] input->stepvol
		; [EAX + 36] input->wscale
		; ...
		MOVSS XMM0, [EAX+28]
		MOVSS [stepind], XMM0
		prints msg
		printss stepind
		prints nl

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante

global addVettori
	
	dim equ 28
	inizio2 equ 24
	inizio1 equ 20
	ris equ 16
	v2 equ 12
	v1 equ 8
	

addVettori:
	
    push		ebp
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push		ebx		; salva i registri da preservare
	push		esi
	push		edi



	
	XOR ESI,ESI
    MOV EDI, [EBP+dim]
    MOV EBX,[EBP+v1]
    MOV ECX,[EBP+v2]
    MOV EDX,[EBP+ris]
    MOV EAX,[EBP+inizio1]
    SHL EAX, 2
    ADD EBX,EAX
    MOV EAX,[EBP+inizio2]
    SHL EAX, 2
    ADD EBX,EAX
	cicloaddvettori:SUB EDI,4
			CMP EDI,0
			JL fineAddVettori
			MOVUPS XMM0,[EBX + 4*ESI]
			ADDPS XMM0,[ECX+4*ESI]
			MOVUPS [EDX+4*ESI],XMM0
			ADD ESI,4	           
        jmp cicloaddvettori
            
	fineAddVettori: ADD EDI,4
                    
	ciclofinevettori:SUB EDI,1
			CMP EDI,0
			JL e1
			MOVSS XMM0,[EBX+4*ESI]
			ADDSS XMM0,[ECX+4*ESI]
			MOVSS [EDX+4*ESI],XMM0
			ADD ESI,1
			JMP ciclofinevettori
   
	e1:
    pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret			; torna alla funzione C chiamante


global subVettori
	
	dim equ 28
	inizio2 equ 24
	inizio1 equ 20
	ris equ 16
	v2 equ 12
	v1 equ 8
	

subVettori:
	
    push		ebp
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push		ebx		; salva i registri da preservare
	push		esi
	push		edi



	
	XOR ESI,ESI
    MOV EDI, [EBP+dim]
    MOV EBX,[EBP+v1]
    MOV ECX,[EBP+v2]
    MOV EDX,[EBP+ris]

    MOV EAX,[EBP+inizio1]
    SHL EAX, 2
    ADD EBX,EAX
    MOV EAX,[EBP+inizio2]
    SHL EAX, 2
    ADD EBX,EAX
	ciclosubvettori:SUB EDI,4
			CMP EDI,0
			JL fineSubVettori
			MOVUPS XMM0,[EBX + 4*ESI]
			SUBPS XMM0,[ECX+4*ESI]
			MOVUPS [EDX+4*ESI],XMM0
			ADD ESI,4	           
        jmp ciclosubvettori
            
	fineSubVettori: ADD EDI,4
                    
	ciclofinesubvettori:SUB EDI,1
			CMP EDI,0
			JL e2
			MOVSS XMM0,[EBX+4*ESI]
			SUBSS XMM0,[ECX+4*ESI]
			MOVSS [EDX+4*ESI],XMM0
			ADD ESI,1
			JMP ciclofinesubvettori
   
	e2:
    pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret			; torna alla funzione C chiamante
		


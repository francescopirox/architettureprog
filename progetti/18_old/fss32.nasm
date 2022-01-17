; ---------------------------------------------------------
; Regressione con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;
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

		
;global addVettori:
	dim equ 8
	inizio2 equ 12
	inizio1 equ 16
	ris equ 20
	v2 equ 24
	v1 equ 28
	

addVettori:
	
	push		ebp		; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push		ebx		; salva i registri da preservare
	push		esi
	push		edi
		
		
	XOR ESI,ESI
	MOV EAX,[EBP+dim]
	cicloaddvettori:SUB EAX,4
			CMP EAX,0
			JL fineAddVettori
			MOVAPS XMM0,[EBP+v1+4*ESI]
			ADDPS XMM0,[EBP+v2+4*ESI]
			MOVAPS [EBP+ris+4*ESI],XMM0
			ADD ESI,4
			jmp cicloaddvettori
	fineAddVettori: ADD EAX,4
	ciclofinevettori:SUB EAX,1
			CMP EAX,0
			JL e1
			MOVSS XMM0,[EBP+v1+inizio1+4*ESI]
			ADDSS XMM0,[EBP+v2+inizio2+4*ESI]
			MOVSS [EBP+ris+4*ESI],XMM0
			ADD ESI,1
			JMP ciclofinevettori
	e1: pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret			; torna alla funzione C chiamante
		
			

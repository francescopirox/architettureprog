; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
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
;     nasm -f elf64 regression64.nasm
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

section .bss			; Sezione contenente dati non inizializzati

alignb 32
stepind		resq		1
alignb 32
de		resq		1
alignb 32
risSommaEu	resq		1

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
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

; ------------------------------------------------------------
; Funzione prova
; ------------------------------------------------------------
global prova

msg	db 'stepind:',0
nl	db 10,0

prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

		; ------------------------------------------------------------
		; I parametri sono passati nei registri
		; ------------------------------------------------------------
		; rdi = indirizzo della struct input

		; esempio: stampa input->stepind
		; RDI contiente l'indirizzo della struttura contenente i parametri
		; [EDI]	input->x
		; [EDI + 8] input->xh
		; [EDI + 16] input->c
		; [EDI + 24] input->r
		; [EDI + 32] input->nx
		; [EDI + 36] input->d
		; [EDI + 40] input->iter
		; [EDI + 48] input->stepind
		; [EDI + 56] input->stepvol
		; [EDI + 64] input->wscale
		; ...
		VMOVSD		XMM0, [RDI+48]
		VMOVSD		[stepind], XMM0
		prints 		msg
		printsd		stepind
		prints 		nl
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		; ripristina il Base Pointer
		ret				; torna alla funzione C chiamante


global distEuclideaAsm 
	
	
distEuclideaAsm:
push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali
		
		;inserire cicli
		
        VXORPD YMM1,YMM1
	XOR RAX,RAX
        XOR RBX,RBX
        SHL RDX,3
        SHL RCX,3
        ADD RDI,RDX
        ADD RSI,RCX
        
cicloSommaEuclidea:	
	SUB R8,4
	CMP R8,0
	JL fineSommaEuclidea	       
        VMOVAPD YMM2,[RSI + RAX]        
        VSUBPD YMM2,[RDI +RAX]
        VMULPD YMM2,YMM2
	VADDPD YMM1,YMM2
       	ADD RAX,32
       	JMP cicloSommaEuclidea

fineSommaEuclidea:
	ADD R8,3
	VXORPD XMM3,XMM3
cicloFineSommaEuclidea:	
	CMP R8,0
	JL e1
        VMOVSD XMM2,[RSI+RAX]
        VSUBSD XMM2,[RDI+RAX]
        VMULSD XMM2,XMM2,XMM2
	VADDSD XMM3,XMM2
	SUB R8,1
	ADD RAX,8
	JMP cicloFineSommaEuclidea
	
e1:     
	VHADDPD YMM1,YMM1
        VEXTRACTF128 XMM0,YMM1,1b ;<-------
        VADDPD XMM0,XMM1
        VADDPD XMM0,XMM3
        VSQRTSD XMM0,XMM0
        VMOVSD [risSommaEu],XMM0 
       
      
	popaq				; ripristina i registri generali
	mov		rsp, rbp	; ripristina lo Stack Pointer
	pop		rbp		; ripristina il Base Pointer
	VMOVSD XMM0,[risSommaEu] 
	ret				; torna alla funzione C chiamante
	

global pesoTot 
	
	
pesoTot:
push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali
		
		;inserire cicli
		
        VXORPD YMM1,YMM1
        XOR RAX,RAX
cicloPesoTot:	
	SUB RSI,4
	CMP RSI,0
	JL finePesoTot	       
        VADDPD YMM1,[RDI + RAX]        
       	ADD RAX,32
       	JMP cicloPesoTot

finePesoTot:
	ADD RSI,3
	VXORPD XMM2,XMM2
	
cicloFinePesoTot:	
	CMP RSI,0
	JL e2
        VADDSD XMM2,[RDI+RAX]	
	SUB RSI,1
	ADD RAX,8
	JMP cicloFinePesoTot
	
e2:     
	VHADDPD YMM1,YMM1
        VEXTRACTF128 XMM0,YMM1,1b ;<-------
        VADDPD XMM0,XMM1
        VADDPD XMM0,XMM2
        VMOVSD [de],XMM0 
       
      
	popaq				; ripristina i registri generali
	mov		rsp, rbp	; ripristina lo Stack Pointer
	pop		rbp		; ripristina il Base Pointer
	VMOVSD XMM0,[de] 
	ret				; torna alla funzione C chiamante

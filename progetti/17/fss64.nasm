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
alignb 32
prova1	resq		4
alignb 32
cv	resq		1

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


global distEuclidea 
	
	
distEuclidea:
push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali
		
		;inserire cicli
		
        VXORPD YMM1,YMM1
	XOR RAX,RAX
        XOR RBX,RBX
        
        
cicloSommaEuclidea:	
	SUB RDX,4
	;CMP RDX,0
	JL fineSommaEuclidea	       
        VMOVAPD YMM2,[RSI + RAX]        
        VSUBPD YMM2,[RDI + RAX]
        VMULPD YMM2,YMM2
	VADDPD YMM1,YMM2
       	ADD RAX,32
       	JMP cicloSommaEuclidea

fineSommaEuclidea:
	ADD RDX,3
	VXORPD XMM3,XMM3
cicloFineSommaEuclidea:	
	CMP RDX,0
	JL e1
        VMOVSD XMM2,[RSI+RAX]
        VSUBSD XMM2,[RDI+RAX]
        VMULSD XMM2,XMM2,XMM2
	VADDSD XMM3,XMM2
	SUB RDX,1
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
	
	
;---------------------------------

global copyAlnVectorAsm

copyAlnVectorAsm:
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq	
	
    	
    	MOV	R8,	RSI
    	AND	R8,	3
    	JZ	modulo4
    	prints msg
    	PUSH 	RDI
    	PUSH 	RSI
    	PUSH 	RDX
    	getmem 8, RDX
    	MOV	[cv],	RAX
    	POP	RDX
    	POP	RSI
    	POP	RDI
    	
    	SHL 	RDX,	3
	XOR	RBX,	RBX
    	SHL	RSI,	3
    	ADD	RDI,	RSI
    	
	
    		
cicloquozientecopy:
	
    	SUB RDX,32					; dim=dim-4
	;CMP EDI,0
	JL finecicloquozientecopy
	VMOVUPD	YMM0,	[RDI+RBX]
	VMOVAPD	[RAX+RBX],	YMM0
	ADD	RBX,	32
	JMP	cicloquozientecopy
	
finecicloquozientecopy:
	ADD	RDX,	32
	ciclorestocopy:
		SUB 	RDX,	8
		JL	ecopy
		VMOVSS	XMM0,	[RDI+RBX]
		VMOVSS	[RAX+RBX],	XMM0
		ADD 	RBX,	8
		JMP	ciclorestocopy
modulo4:
	
    	SHL	RSI,	3
    	ADD	RDI,	RSI
	MOV	[cv],	RDI	
	
ecopy:
	 
	popaq				; ripristina i registri generali
	mov		rsp, rbp	; ripristina lo Stack Pointer
	pop		rbp		; ripristina il Base Pointer
	MOV 	RAX,	[cv]
	ret				; torna alla funzione C chiamante



;--------------------
global addVettori

addVettori:
	push		rbp				
	mov			rbp, rsp			
	pushaq	
	
	;RDI=v1
	;RSI=v2
	;RDX=ris
	;RCX=dim
	
	XOR 			RAX,RAX			;i=0
	
cicloaddVettori:
	SUB 			RCX, 4			;dim=dim-4
	CMP 			RCX, 0			;dim==0?
	JL 			fineaddVettori
	VMOVAPD 	YMM1, [RDI + RAX] 	;YMM1=[V1, V1, V1, V1]
	VADDPD 		YMM1, [RSI + RAX]	;YMM1=[V1+V2, V1+V2, V1+V2, V1+V2]
	VMOVAPD 	[RDX + RAX], YMM1	; metto nel vettore ris i 4 double
	ADD 		RAX, 32			; i+=4
	JMP 			cicloaddVettori
	
fineaddVettori:
	ADD 		RCX, 3			;dim +=4? o 3?
	;VXORPD 		YMM1,YMM1		;YMM1=[0, 0, 0, 0]
	
ciclofineaddVettori:
	CMP 			RCX, 0			;dim==0?
	JL 			e3
	VMOVSD 		XMM1, [RDI + RAX] 	;YMM1=[0, 0, 0, V1]
	VADDSD 		XMM1, [RSI + RAX]	;YMM1=[0, 0, 0, V1+V2]
	VMOVSD 		[RDX + RAX], XMM1  ;metto il risultato in ris
	SUB 			RCX, 1			;dim--
	ADD 		RAX, 8			;i++
	JMP 			ciclofineaddVettori
	
e3:
	popaq					
	mov			rsp, rbp	
	pop			rbp
	ret
;---------------------------------------------------
global subVettori

subVettori:
	push		rbp				
	mov			rbp, rsp			
	pushaq	
	
	;RDI=v1
	;RSI=v2
	;RDX=ris
	;RCX=dim
	
	XOR 			RAX,RAX			;i=0
	
ciclosubVettori:
	SUB 			RCX, 4			;dim=dim-4
	CMP 			RCX, 0			;dim==0?
	JL 			finesubVettori
	VMOVAPD 	YMM1, [RDI + RAX] 	;YMM1=[V1, V1, V1, V1]
	VSUBPD 		YMM1, [RSI + RAX]	;YMM1=[V1-V2, V1-V2, V1-V2, V1-V2]
	VMOVAPD 	[RDX + RAX], YMM1	; metto nel vettore ris i 4 double
	ADD 		RAX, 32			; i+=4
	JMP 			ciclosubVettori
	
finesubVettori:
	ADD 		RCX, 3			;dim +=4? o 3?
	;VXORPD 		YMM1,YMM1		;YMM1=[0, 0, 0, 0]
	
ciclofinesubVettori:
	CMP 			RCX, 0			;dim==0?
	JL 			e4
	VMOVSD 		XMM1, [RDI + RAX] 	;YMM1=[0, 0, 0, V1]
	VSUBSD 		XMM1, [RSI + RAX]	;YMM1=[0, 0, 0, V1-V2]
	VMOVSD 		[RDX + RAX], XMM1  ;metto il risultato in ris
	SUB 			RCX, 1			;dim--
	ADD 		RAX, 8			;i++
	JMP 			ciclofinesubVettori
	
e4:
	popaq					
	mov			rsp, rbp	
	pop			rbp
	ret

;------------------------------------------------------------------

global prodVet_x_Scalare

prodVet_x_Scalare:
	push		rbp
	mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
	pushaq							; salva i registri generali
	
	; RDI = v1
	; XMM0 = S
	; RSI = ris
	; RDX = dim
	
	XOR 		RAX,RAX
	;VBROADCASTSD	YMM0, XMM0
	
	;VSHUFPD 	YMM0, XMM0, 00000000b	; XMM0 = [S, S, S, S]
	VXORPD 	YMM3, YMM3
	
cicloProdVet_x_Scalare:
	SUB 		RDX,4					; DIM-4
	CMP 		RDX,0					; DIM == 0?
	JL 		fineProdVet_x_Scalare		; se si jumpa all'ultima iterazione
	VMOVAPD YMM3, [RDI + RAX]			; YMM3 = [V1, V1, V1, V1]
	VMULPD	YMM3, YMM0				; YMM3 = [V1*S, V1*S, V1*S, V1*S]
	VMOVAPD [RSI+RAX], YMM3			; LO RIMETTO A POSTO in memoria
	ADD 	RAX,32					; i= i+4
	JMP 		cicloProdVet_x_Scalare

fineProdVet_x_Scalare:
	ADD 	RDX,3					; dim=dim+4
	VXORPD 	YMM3,YMM3
	
ciclofineProdVet_x_Scalare:
	CMP 		RDX,0					; dim == 0 ?
	JL 		e5						; se si fine
    VMOVSD 	XMM3, [RDI + RAX]			; YMM3 = [0, 0, 0, V1]
	VMULSD	XMM3, XMM0				; YMM3 = [0*S, 0*S, 0*S, V1*S]
	VMOVSD 	[RSI+RAX], XMM3			; LO RIMETTO A POSTO in memoria
	SUB 		RDX,1					; dim--
	ADD 	RAX,8					; i++
	JMP 		ciclofineProdVet_x_Scalare

e5:
	popaq							; ripristina i registri generali
	mov		rsp, rbp					; ripristina lo Stack Pointer
	pop		rbp						; ripristina il Base Pointer
			
	ret								; torna alla funzione C chiamante
;--------------------------------------------------

global prodScalare

prodScalare:
	push		rbp
	mov			rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
	pushaq					; salva i registri generali
	
	; XMM0 prima ora RDI  = v1
	; XMM1 prima ora RSI  = v2
	; RDX = dim
	
	XOR 		RAX,RAX
	VXORPD 	YMM1, YMM1
	VXORPD 	YMM2, YMM2			;contatore
	
cicloProdScalare:
	SUB 		RDX,4						; DIM-4
	CMP 		RDX,0						; DIM == 0?
	JL 		fineProdScalare				; se si jumpa all'ultima iterazione
	VMOVAPD YMM0, [RDI + RAX]				; YMM0 = [V1, V1, V1, V1]
	VMULPD	YMM0, [RSI + RAX]				; YMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
	VADDPD	YMM1, YMM0
	ADD 	RAX,32						; i= i+4
	JMP 		cicloProdScalare

fineProdScalare:
	ADD 	RDX,3						; dim=dim+4
	VXORPD 	YMM0, YMM0
	
ciclofineProdScalare:
	CMP 		RDX,0						; dim == 0 ?
	JL 		e6							; se si fine
        VMOVSD 	XMM0, [RDI + RAX]				; YMM0 = [0, 0, 0, V1]
	VMULSD	XMM0, [RSI + RAX]				; YMM0 = [0*V2, 0*V2, 0*V2, V1*V2]
	VADDSD	XMM2, XMM0				
	SUB 		RDX,1						; dim--
	ADD 	RAX,8						; i++
	JMP 		ciclofineProdScalare

e6:
	VHADDPD 	YMM1, YMM1
	VEXTRACTF128 	XMM3, YMM1, 1b 
	
	VADDSD 		XMM1, XMM3
	VADDSD 		XMM2, XMM1
	;VEXTRACTF128 XMM0, YMM2, 1b 			;<-------
	
	VMOVSD 		[risSommaEu], XMM2

	popaq								; ripristina i registri generali
	mov			rsp, 	rbp					; ripristina lo Stack Pointer
	pop			rbp						; ripristina il Base Pointer
	VMOVSD 		XMM0, [risSommaEu] 
	ret									; torna alla funzione C chiamante


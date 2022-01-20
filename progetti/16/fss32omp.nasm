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
	alignb 16
	retps		resd		1
	alignb 16
	ss1		resd		4
	alignb 16
	retpt		resd		1


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
	
	dim equ 20
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

    cicloaddvettori8:SUB EDI,32
			;CMP EDI,0
			JL cicloaddvettoriold
			MOVAPS XMM0,[EBX + 4*ESI]
			ADDPS XMM0,[ECX+4*ESI]
			MOVAPS [EDX+4*ESI],XMM0
           
            MOVAPS XMM1,[EBX + 4*ESI +16]
			ADDPS XMM1,[ECX+4*ESI +16]
			MOVAPS [EDX+4*ESI+16],XMM1
			
            MOVAPS XMM0,[EBX + 4*ESI+32]
			ADDPS XMM0,[ECX+4*ESI+32]
			MOVAPS [EDX+4*ESI+32],XMM0
           
            MOVAPS XMM1,[EBX + 4*ESI+48]
			ADDPS XMM1,[ECX+4*ESI+48]
			MOVAPS [EDX+4*ESI+48],XMM1
			
            MOVAPS XMM0,[EBX + 4*ESI+64]
			ADDPS XMM0,[ECX+4*ESI+64]
			MOVAPS [EDX+4*ESI+64],XMM0
            
            MOVAPS XMM1,[EBX + 4*ESI+80]
			ADDPS XMM1,[ECX+4*ESI+80]
			MOVAPS [EDX+4*ESI+80],XMM1
			
            MOVAPS XMM0,[EBX + 4*ESI+96]
			ADDPS XMM0,[ECX+4*ESI+96]
			MOVAPS [EDX+4*ESI+96],XMM0
            
            MOVAPS XMM1,[EBX + 4*ESI+112]
			ADDPS XMM1,[ECX+4*ESI+112]
			MOVAPS [EDX+4*ESI+112],XMM1
			ADD ESI,32	 	      	           
        jmp cicloaddvettori8
    
    cicloaddvettoriold:
            ADD EDI, 32
	cicloaddvettori2:SUB EDI,8
			;CMP EDI,0
			JL fineAddVettori
			MOVAPS XMM0,[EBX + 4*ESI]
			ADDPS XMM0,[ECX+4*ESI]
			MOVAPS [EDX+4*ESI],XMM0
            ;ADD ESI,4
            MOVAPS XMM1,[EBX + 4*ESI+16]
			ADDPS XMM1,[ECX+4*ESI+16]
			MOVAPS [EDX+4*ESI+16],XMM1
			ADD ESI,8	           
        jmp cicloaddvettori2
            
	fineAddVettori: ADD EDI,8
                    
	ciclofinevettori:SUB EDI,1
			;CMP EDI,0
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
	
	dim equ 20
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
    MOV EDI,[EBP+dim]
    MOV EBX,[EBP+v1]
    MOV ECX,[EBP+v2]
    MOV EDX,[EBP+ris]
    
    	ciclosubvettori8:SUB EDI,32
			;CMP EDI,0
			JL ciclosubvettoriold
			MOVAPS XMM0,[EBX + 4*ESI]
			SUBPS XMM0,[ECX+4*ESI]
            ;SUBPS XMM0,XMM1
			MOVAPS [EDX+4*ESI],XMM0
           
            MOVAPS XMM1,[EBX + 4*ESI +16]
			SUBPS XMM1,[ECX+4*ESI +16]
            ;SUBPS XMM0,XMM1
			MOVAPS [EDX+4*ESI +16],XMM1
		
            MOVAPS XMM0,[EBX + 4*ESI +32]
			SUBPS XMM0,[ECX+4*ESI+32]
            ;SUBPS XMM0,XMM1
			MOVAPS [EDX+4*ESI +32],XMM0
           
            MOVAPS XMM1,[EBX + 4*ESI+48]
			SUBPS XMM1,[ECX+4*ESI+48]
            ;SUBPS XMM0,XMM1
			MOVAPS [EDX+4*ESI+48],XMM1
			
            MOVAPS XMM0,[EBX + 4*ESI+64]
			SUBPS XMM0,[ECX+4*ESI+64]
            ;SUBPS XMM0,XMM1
			MOVAPS [EDX+4*ESI+64],XMM0
           
            MOVAPS XMM1,[EBX + 4*ESI+80]
			SUBPS XMM1,[ECX+4*ESI+80]
            ;SUBPS XMM0,XMM1
			MOVAPS [EDX+4*ESI+80],XMM1
			
            MOVAPS XMM0,[EBX + 4*ESI+96]
			SUBPS XMM0,[ECX+4*ESI+96]
            ;SUBPS XMM0,XMM1
			MOVAPS [EDX+4*ESI+96],XMM0
          
            MOVAPS XMM1,[EBX + 4*ESI+112]
			SUBPS XMM1,[ECX+4*ESI+112]
            ;SUBPS XMM0,XMM1
			MOVAPS [EDX+4*ESI+112],XMM1
			ADD ESI,32	   	     	   	           
        jmp ciclosubvettori8
    ciclosubvettoriold:
            ADD EDI,32
	ciclosubvettori2:SUB EDI,8
			;CMP EDI,0
			JL fineSubVettori
			MOVAPS XMM0,[EBX + 4*ESI]
			SUBPS XMM0,[ECX+4*ESI]
            ;SUBPS XMM0,XMM1
			MOVAPS [EDX+4*ESI],XMM0
            ;ADD ESI,4
            MOVAPS XMM1,[EBX + 4*ESI+16]
			SUBPS XMM1,[ECX+4*ESI+16]
            ;SUBPS XMM0,XMM1
			MOVAPS [EDX+4*ESI+16],XMM1
			ADD ESI,8	           
        jmp ciclosubvettori2
            
	fineSubVettori: ADD EDI,8
                    
	ciclofinesubvettori:SUB EDI,1
			;CMP EDI,0
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



;--------------------------

global prodVet_x_Scalare

	dimV equ 20
	ris equ 16
	s equ 12
	v equ 8

prodVet_x_Scalare:
	push		ebp
	mov		ebp, esp		; il Base Pointer punta al Record di Attivazione corrente
	push		ebx		; salva i registri da preservare
	push		esi
	push		edi

	XOR		ESI,  ESI					; i=0
	MOV		EDI,	[EBP+dimV]			; EDI = dim
	MOV 	EBX, [EBP+ris]			; EBX = RIS (� L'indirizzo di UN VETTORE)
	MOVSS 	XMM0, [EBP+s]			; XMM0 = S  [0,0,0,S] (� UNO SCALARE)
	MOV 	EDX, [EBP+v]			; EDX =V1 (VETTORE DI PARTENZA)
	
	;XORPS	XMM0, XMM0
	SHUFPS 	XMM0, XMM0, 00000000b	; XMM0 = [S, S, S, S]

cicloProdVet_x_Scalare8:
	SUB 		EDI,32					; DIM-4
	;CMP 		EDI,0					; DIM == 0?
	JL 		cicloProdVet_x_Scalareold		; se si jumpa all'ultima iterazione
	
    MOVAPS 	XMM1,[EDX + 4*ESI]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM1, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI], XMM1
	

    MOVAPS 	XMM2,[EDX + 4*ESI+16]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM2, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI+16], XMM2
	

    MOVAPS 	XMM3,[EDX + 4*ESI+32]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM3, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI+32], XMM3
	

    MOVAPS 	XMM4,[EDX + 4*ESI+48]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM4, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI+48], XMM4
				

    MOVAPS 	XMM5,[EDX + 4*ESI +64]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM5, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI+ 64], XMM5
					; i= i+4

    MOVAPS 	XMM6,[EDX + 4*ESI + 80]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM6, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI + 80], XMM6


    MOVAPS 	XMM5,[EDX + 4*ESI +96]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM5, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI + 96], XMM5


    MOVAPS 	XMM6,[EDX + 4*ESI +112]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM6, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI +112], XMM6
	ADD 	ESI,32
					
	JMP 	cicloProdVet_x_Scalare8
cicloProdVet_x_Scalareold:
    ADD EDI,    32
cicloProdVet_x_Scalare2:
	SUB 		EDI,8					; DIM-4
	;CMP 		EDI,0					; DIM == 0?
	JL 		fineProdVet_x_Scalare		; se si jumpa all'ultima iterazione
	MOVAPS 	XMM1,[EDX + 4*ESI]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM1, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI], XMM1

    MOVAPS 	XMM2,[EDX + 4*ESI+16]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM2, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVUPS 	[EBX+ 4*ESI+16], XMM2
	ADD 	ESI,8					; i= i+4
	JMP 	cicloProdVet_x_Scalare2

fineProdVet_x_Scalare:
	ADD 	EDI,8					; dim=dim+4

ciclofineProdVet_x_Scalare:
	SUB 		EDI,1					; dim=dim-1
	;CMP 		EDI,0					; dim==0?
	JL 		e3						; se si ho finito e faccio le pop dei registri dallo stack
	;XORPS	XMM1, XMM1				; XMM1 = [0,0,0,0]
	MOVSS 	XMM1,[EDX+4*ESI]			; XMM1 = [0,0,0,v1[i]]
	MULSS 	XMM1, XMM0				; XMM1 = [0*S, 0*S, 0*S, v1[i]*S]
	MOVSS 	[EBX+4*ESI], XMM1
	INC 		ESI						; i++
	JMP 		ciclofineProdVet_x_Scalare

e3:
	pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret			; torna alla funzione C chiamante
;--------------------------------------------------

global prodScalare
	
	dimps equ 16
	v2 equ 12
	v1 equ 8
	
prodScalare:
	
	push		ebp
	mov		ebp, esp		; il Base Pointer punta al Record di Attivazione corrente
	push		ebx		; salva i registri da preservare
	push		esi
	push		edi
	
	XOR	ESI, 	ESI			; i=0
	MOV	EDI,	[EBP+dimps]		; EDI = dim
	MOV 	EBX,	[EBP+v1]		; EBX=V1
	MOV 	ECX,	[EBP+v2]		; ECX=V2
; EDX=RIS  (� un float)

	
	XORPS 	XMM1,	XMM1			; XMM1=0 IL MIO CONTATORE DOVE OGNI V

cicloprodScalare8:	
	SUB 	EDI,    16				; DIM-4
	;CMP 	EDI,    0				; DIM == 0?
	JL 		cicloprodScalarold		; se si jumpa all'ultima iterazione
	MOVAPS 	XMM0,	[EBX + 4*ESI]		; XMM0 = [V1, V1, V1, V1]
	MULPS	XMM0,	[ECX + 4*ESI]		; XMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
	ADDPS 	XMM1,	XMM0					;
    
    MOVAPS 	XMM2,	[EBX + 4*ESI+16]		; XMM0 = [V1, V1, V1, V1]
	MULPS	XMM2,	[ECX + 4*ESI+16]		; XMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
	ADDPS 	XMM1,	XMM2
    
    MOVAPS 	XMM3,	[EBX + 4*ESI +32]		; XMM0 = [V1, V1, V1, V1]
	MULPS	XMM3,	[ECX + 4*ESI +32]		; XMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
	ADDPS 	XMM1,	XMM3					;i= i+4

    MOVAPS 	XMM4,	[EBX + 4*ESI+48]		; XMM0 = [V1, V1, V1, V1]
	MULPS	XMM4,	[ECX + 4*ESI+48]		; XMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
	ADDPS 	XMM1,	XMM4

	ADD 	ESI,16					;i= i+4
    JMP cicloprodScalare8

cicloprodScalarold:
	ADD EDI, 16
cicloprodScalare2:	
	SUB 	EDI,    8				; DIM-4
	;CMP 	EDI,    0				; DIM == 0?
	JL 		fineProdScalare			; se si jumpa all'ultima iterazione
	MOVAPS 	XMM0,	[EBX + 4*ESI]		; XMM0 = [V1, V1, V1, V1]
	MULPS	XMM0,	[ECX + 4*ESI]		; XMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
	ADDPS 	XMM1,	XMM0					;i= i+4
    MOVAPS 	XMM2,	[EBX + 4*ESI+16]		; XMM0 = [V1, V1, V1, V1]
	MULPS	XMM2,	[ECX + 4*ESI+16]		; XMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
	ADDPS 	XMM1,	XMM2
	ADD 	ESI,8					;i= i+4
    JMP cicloprodScalare2
            
fineProdScalare: 
	ADD 	EDI,	8					;dim=dim+4
	XORPS	XMM3,	XMM3    
ciclofineProdScalare:	
	SUB 	EDI,	1					;dim=dim-1
	;CMP 	EDI,	0					;dim==0?
	JL 	e4						;se si ho finito e faccio le pop dei registri dallo stack 
	MOVSS 	XMM0,	[EBX+4*ESI]			;XMM0 = [0, 0, 0,	v1[i]]
	MULSS 	XMM0,	[ECX+4*ESI]			;XMM0 = [0, 0, 0, v1[i]*v2[i]]
	ADDSS	XMM3,	XMM0				;XMM1 = XMM1 + v1[i]*v2[i]
	INC 		ESI						;i++
	JMP 		ciclofineProdScalare
   
e4:	
	
	HADDPS 	XMM1,	XMM1
	HADDPS 	XMM1,	XMM1
	
    	ADDSS	XMM1,	XMM3				
	MOVSS	[retps],	XMM1
    	FLD	dword [retps]
	pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret			; torna alla funzione C chiamante
	
;----------------------------------------------------------------------------
global pesoTot

       dimp equ 12
       v equ 8

pesoTot: 
	push   ebp
    	mov    ebp, esp
    	push   ebx
    	push   esi
    	push   edi

	XOR ESI,ESI					; i=0
    	MOV EDI, [EBP+dimp]			; EDI = dim
    	MOV EBX, [EBP+v]				; EBX = v
    	XORPS XMM0, XMM0				; XMM0 = [0,0,0,0]
    	XORPS XMM1, XMM1				; XMM1 = [0,0,0,0]

cicloPesoTot4: 
	SUB EDI,16					; dim=dim-4
	;CMP EDI,0					; dim==0?
    JL  	cicloPesoTotOld					; se si vado a stop
	ADDPS XMM0, [EBX+ ESI*4]		; XMM0 = [v,v,v,v]
    					; i=i+4
    ADDPS XMM0, [EBX+ ESI*4+16]		; XMM0 = [v,v,v,v]
    ADDPS XMM0, [EBX+ ESI*4+32]		; XMM0 = [v,v,v,v]
    					; i=i+4
    
    ADDPS XMM0, [EBX+ ESI*4+112]		; XMM0 = [v,v,v,v]
    ADD ESI, 16					; i=i+4
	JMP cicloPesoTot4	
	

cicloPesoTotOld:
    ADD EDI,16
cicloPesoTot: 
	SUB EDI,8					; dim=dim-4
	;CMP EDI,0					; dim==0?
    JL  	finePesoTot					; se si vado a stop
	ADDPS XMM0, [EBX+ ESI*4]		; XMM0 = [v,v,v,v]
    					; i=i+4
    ADDPS XMM0, [EBX+ ESI*4+16]		; XMM0 = [v,v,v,v]
    ADD ESI, 8					; i=i+4
	JMP cicloPesoTot	
	
finePesoTot: 
	ADD EDI, 8
	
ciclofinePesoTot: 
	SUB EDI,1					; dim--
	CMP EDI, 0					; dim ==0?
	JL e5						; se si ho finito
      	ADDSS XMM1, [EBX+ESI*4]		; XMM1 = [0,0,0,v]
	INC ESI						; i++
     	JMP ciclofinePesoTot
e5:
	HADDPS XMM0, XMM0		
	HADDPS XMM0,XMM0			; prima somma parziale (multipli di 4)
	ADDSS XMM0, XMM1			; somma delle due somme parziali

    	MOVSS [retpt], XMM0

	FLD dword [retpt]
   	pop edi 
    	pop esi
    	pop ebx
    	mov esp, ebp 
    	pop ebp 
    	ret
;-------------------------
global copyAlnVector

vcopy	equ	8
inizio	equ	12
dimcopy	equ	16

;msg	db	'stepind:',32,0

copyAlnVectorAsm:
	push   ebp
    	mov    ebp, esp
    	push   ebx
    	push   esi
    	push   edi
    	
    	MOV	EDX,	[EBP+inizio]
    	AND	EDX,	3
    	JZ	modulo4
    	
    	prints msg
    	getmem 4,[EBP+dim]
    	
	XOR	ESI,	ESI				; i=0
    	MOV	EDI,	[EBP+dimcopy]			; EDI = dim
    	MOV	ECX,	[EBP+inizio]
    	SHL	ECX,	2
    	MOV	EBX,	[EBP+vcopy]
    	ADD	EBX,	ECX
    	
	
    		
cicloquozientecopy:
	
    	SUB EDI,4					; dim=dim-4
	;CMP EDI,0
	JL finecicloquozientecopy
	;MOVUPS	XMM0,	[EBX+ESI]
	;MOVAPS	[EAX+ESI],	XMM0
	ADD	ESI,	16
	JMP	cicloquozientecopy
	
finecicloquozientecopy:
	ADD	EDI,	4
	ciclorestocopy:
		SUB 	EDI,	1
		;CMP	EDI,	0
		JL	ecopy
		MOVSS	XMM0,	[ESI+EAX]
		MOVSS	[ESI+EBX],	XMM0
		ADD 	ESI,	4
		JMP	ciclorestocopy
modulo4:
    ;prints msg
	MOV	ECX,	[EBP+inizio]
    SHL	ECX,	2
    MOV	EBX,	[EBP+vcopy]
    ADD	EBX,	ECX
	MOV	EAX,	EBX	
	
ecopy:
	pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret	


;------------------
global distEuclidea
	dimeucl equ 16
       	v2 equ 12
       	v1 equ 8

distEuclidea: 
	    push   ebp
    	mov    ebp, esp
    	push   ebx
    	push   esi
    	push   edi

	XOR ESI,ESI
    	MOV EDI, [EBP+dimeucl]
	MOV EAX, [EBP+v1]
    	MOV EBX, [EBP+v2]
    	XORPS XMM1, XMM1  ;registro per le somme a multipli di 4
cicloDistEuclidea4: 
	SUB EDI,32
    CMP EDI,0
    JL  cicloDistEuclideaOld

    MOVAPS XMM0,[EAX+ESI*4]      ;XMM0=[v1, v1, v1, v1]
	SUBPS XMM0, [EBX+ESI*4]        ;XMM0=[v1-v2, v1-v2, v1-v2, v1-v2]
	MULPS XMM0, XMM0  ;XMM0=[v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2]
	ADDPS XMM1, XMM0    ;sommo le differenze dei quadrati sulle celle di XMM1
   
    MOVAPS XMM2,[EAX+ESI*4+16]      ;XMM0=[v1, v1, v1, v1]
    SUBPS XMM2, [EBX+ESI*4+16]        ;XMM0=[v1-v2, v1-v2, v1-v2, v1-v2]
	MULPS XMM2, XMM2  ;XMM0=[v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2]
	ADDPS XMM1, XMM2    ;sommo le differenze dei quadrati sulle celle di XMM1
    
    MOVAPS XMM0,[EAX+ESI*4 +32]      ;XMM0=[v1, v1, v1, v1]
	SUBPS XMM0, [EBX+ESI*4 +32]        ;XMM0=[v1-v2, v1-v2, v1-v2, v1-v2]
	MULPS XMM0, XMM0  ;XMM0=[v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2]
	ADDPS XMM1, XMM0    ;sommo le differenze dei quadrati sulle celle di XMM1
   
    MOVAPS XMM2,[EAX+ESI*4+16 +32]      ;XMM0=[v1, v1, v1, v1]
    SUBPS XMM2, [EBX+ESI*4+16 +32]        ;XMM0=[v1-v2, v1-v2, v1-v2, v1-v2]
	MULPS XMM2, XMM2  ;XMM0=[v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2]
	ADDPS XMM1, XMM2    ;sommo le differenze dei quadrati sulle celle di XMM1
    
    MOVAPS XMM0,[EAX+ESI*4 + 48]      ;XMM0=[v1, v1, v1, v1]
	SUBPS XMM0, [EBX+ESI*4 + 48]        ;XMM0=[v1-v2, v1-v2, v1-v2, v1-v2]
	MULPS XMM0, XMM0  ;XMM0=[v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2]
	ADDPS XMM1, XMM0    ;sommo le differenze dei quadrati sulle celle di XMM1
   
    MOVAPS XMM2,[EAX+ESI*4+ 64]      ;XMM0=[v1, v1, v1, v1]
    SUBPS XMM2, [EBX+ESI*4+ 64]        ;XMM0=[v1-v2, v1-v2, v1-v2, v1-v2]
	MULPS XMM2, XMM2  ;XMM0=[v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2]
	ADDPS XMM1, XMM2    ;sommo le differenze dei quadrati sulle celle di XMM1
    
    MOVAPS XMM0,[EAX+ESI*4 +96]      ;XMM0=[v1, v1, v1, v1]
	SUBPS XMM0, [EBX+ESI*4 +96]        ;XMM0=[v1-v2, v1-v2, v1-v2, v1-v2]
	MULPS XMM0, XMM0  ;XMM0=[v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2]
	ADDPS XMM1, XMM0    ;sommo le differenze dei quadrati sulle celle di XMM1
   
    MOVAPS XMM2,[EAX+ESI*4+16 +112]      ;XMM0=[v1, v1, v1, v1]
    SUBPS XMM2, [EBX+ESI*4+16 +112]        ;XMM0=[v1-v2, v1-v2, v1-v2, v1-v2]
	MULPS XMM2, XMM2  ;XMM0=[v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2]
	ADDPS XMM1, XMM2    ;sommo le differenze dei quadrati sulle celle di XMM1
    ADD ESI, 32
    JMP cicloDistEuclidea4
cicloDistEuclideaOld:
    ADD EDI,32

cicloDistEuclidea: 
	SUB EDI,8
    CMP EDI,0
    JL  fineDistEuclidea

    MOVAPS XMM0,[EAX+ESI*4]      ;XMM0=[v1, v1, v1, v1]
	SUBPS XMM0, [EBX+ESI*4]        ;XMM0=[v1-v2, v1-v2, v1-v2, v1-v2]
	MULPS XMM0, XMM0  ;XMM0=[v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2]
	ADDPS XMM1, XMM0    ;sommo le differenze dei quadrati sulle celle di XMM1
   
    MOVAPS XMM2,[EAX+ESI*4+16]      ;XMM0=[v1, v1, v1, v1]
    SUBPS XMM2, [EBX+ESI*4+16]        ;XMM0=[v1-v2, v1-v2, v1-v2, v1-v2]
	MULPS XMM2, XMM2  ;XMM0=[v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2]
	ADDPS XMM1, XMM2    ;sommo le differenze dei quadrati sulle celle di XMM1
    ADD ESI, 8
    JMP cicloDistEuclidea
fineDistEuclidea:
    XORPS XMM2, XMM2  ;registro per le somme restanti 
	ADD EDI, 8

ciclofineDistEuclidea: 
	SUB EDI,1
        	CMP EDI, 0
         	JL e7
    MOVSS XMM0, [EAX+ESI*4]  ;XMM0=[0, 0, 0, v1]
	SUBSS XMM0, [EBX+ESI*4]      ;XMM0=[0, 0, 0, (v1-v2)]
	MULSS XMM0, XMM0              ;XMM0=[0, 0, 0, (v1-v2)*(v1-v2)]
	ADDSS XMM2, XMM0              ;XMM2=[0, 0, 0, (v1-v2)*(v1-v2)]
    INC ESI
    JMP ciclofineDistEuclidea
e7:
	HADDPS XMM1, XMM1
	HADDPS XMM1,XMM1   ;somma di tutte le differenze dei quadrati dei multipli di 4
	ADDSS XMM1, XMM0    ;somma delle differenze dei quadrati rimasti
    ADDSS XMM1, XMM2    
   	SQRTSS XMM1, XMM1 

    MOVSS [retps], XMM1

	FLD dword [retps]
   	pop edi 
    	pop esi
    	pop ebx
    	mov esp, ebp 
    	pop  ebp 
    	ret
;_________________        




;--------------------------

global prodVet_x_ScalareUn

	dimV equ 20
	ris equ 16
	s equ 12
	v equ 8

prodVet_x_ScalareUn:
	push		ebp
	mov		ebp, esp		; il Base Pointer punta al Record di Attivazione corrente
	push		ebx		; salva i registri da preservare
	push		esi
	push		edi

	XOR		ESI,  ESI					; i=0
	MOV		EDI,	[EBP+dimV]			; EDI = dim
	MOV 	EBX, [EBP+ris]			; EBX = RIS (� L'indirizzo di UN VETTORE)
	MOVSS 	XMM0, [EBP+s]			; XMM0 = S  [0,0,0,S] (� UNO SCALARE)
	MOV 	EDX, [EBP+v]			; EDX =V1 (VETTORE DI PARTENZA)
	
	;XORPS	XMM0, XMM0
	SHUFPS 	XMM0, XMM0, 00000000b	; XMM0 = [S, S, S, S]

cicloProdVet_x_ScalareUn8:
	SUB 		EDI,32					; DIM-4
	;CMP 		EDI,0					; DIM == 0?
	JL 		cicloProdVet_x_ScalareUnold		; se si jumpa all'ultima iterazione
	
    MOVUPS 	XMM1,[EDX + 4*ESI]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM1, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI], XMM1
	

    MOVUPS 	XMM2,[EDX + 4*ESI+16]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM2, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI+16], XMM2
	

    MOVUPS 	XMM3,[EDX + 4*ESI+32]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM3, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI+32], XMM3
	

    MOVUPS 	XMM4,[EDX + 4*ESI+48]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM4, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI+48], XMM4
				

    MOVUPS 	XMM5,[EDX + 4*ESI +64]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM5, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI+ 64], XMM5
					; i= i+4

    MOVUPS 	XMM6,[EDX + 4*ESI + 80]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM6, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI + 80], XMM6


    MOVUPS 	XMM5,[EDX + 4*ESI +96]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM5, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI + 96], XMM5


    MOVUPS 	XMM6,[EDX + 4*ESI +112]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM6, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI +112], XMM6
	ADD 	ESI,32				
	JMP 	cicloProdVet_x_ScalareUn8
cicloProdVet_x_ScalareUnold:
    ADD EDI,    32
cicloProdVet_x_ScalareUn2:
	SUB 		EDI,8					; DIM-4
	;CMP 		EDI,0					; DIM == 0?
	JL 		fineProdVet_x_ScalareUn		; se si jumpa all'ultima iterazione
	MOVUPS 	XMM1,[EDX + 4*ESI]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM1, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI], XMM1

    MOVUPS 	XMM2,[EDX + 4*ESI+16]		; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM2, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI+16], XMM2
	ADD 	ESI,8					; i= i+4
	JMP 	cicloProdVet_x_ScalareUn2

fineProdVet_x_ScalareUn:
	ADD 	EDI,8					; dim=dim+4

ciclofineProdVet_x_ScalareUn:
	SUB 		EDI,1					; dim=dim-1
	;CMP 		EDI,0					; dim==0?
	JL 		e3u						; se si ho finito e faccio le pop dei registri dallo stack
	;XORPS	XMM1, XMM1				; XMM1 = [0,0,0,0]
	MOVSS 	XMM1,[EDX+4*ESI]			; XMM1 = [0,0,0,v1[i]]
	MULSS 	XMM1, XMM0				; XMM1 = [0*S, 0*S, 0*S, v1[i]*S]
	MOVSS 	[EBX+4*ESI], XMM1
	INC 		ESI						; i++
	JMP 		ciclofineProdVet_x_ScalareUn

e3u:
	pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret			; torna alla funzione C chiamante
;--------------------------------------------------


%pragma win32 prefix _

section .data
    c2 dq 2.0
    c5 dq 5.0
    c3 dq 3.0

section .text
    global f1, f2, f3, f1p, f2p, f3p

;  f1 = 2x + 1 
f1:
    push ebp
    mov ebp, esp
    fld qword[ebp + 8]       ; st0 = x
    fld st0                  ; st0 = x,   st1 = x
    frndint                  ; st0 = (целая часть),   st1 = x
    fxch                     ; st0 = x,   st1 = n
    fsub st0, st1       ; st0 = (остаток) = x - (целая часть), st1=(целая часть)
    f2xm1               ; st0 = 2^(остаток) - 1,   st1 = (целая часть)
    fld1
    faddp               ; st0 = 2^(остаток),  st1 = (целая часть)
    fscale              ; st0 = 2^(остаток) * 2^(целая часть) = 2^x, st1 = n
    fstp st1            ; st0 = 2^x
    fld1
    faddp               ; st0 = 2^x + 1
    pop ebp
    ret

;  f2 = x^5
f2:
    push ebp
    mov ebp, esp
    fld qword[ebp + 8]      ; st0 = x
    fld st0                  ; st0 = x, st1 = x
    fmul st0, st0             ; st0 = x^2, st1 = x
    fmul st0, st0             ; st0 = x^4, st1 = x
    fmulp                     ; st0 = x^5
    pop ebp
    ret

; f3  = (1 - x)/3
f3:
    push ebp
    mov ebp, esp
    fld1                      ; st0 = 1
    fld qword[ebp + 8]      ; st0 = x, st1 = 1
    fsubp                     ; st0 = 1 - x
    fdiv qword[c3]           ; st0 = (1 - x)/3
    pop ebp
    ret

; f1p = 2
f1p:
    push ebp
    mov ebp, esp
    fld qword[ebp + 8]     ; st0 = x
    fld st0
    frndint
    fxch
    fsub st0, st1
    f2xm1
    fld1
    faddp
    fscale              ; st0 = 2^x, st1 = ост
    fstp st1            ; st0 = 2^x
    fldln2           ; st0 = ln(2), st1 = 2^x
    fmulp            ; st0 = 2^x · ln(2)
    pop ebp
    ret

; f2p = 5x^4
f2p:
    push ebp
    mov ebp, esp
    fld qword[ebp + 8]      ; st0 = x
    fmul st0, st0             ; st0 = x^2
    fmul st0, st0             ; st0 = x^4
    fmul qword[c5]           ; st0 = 5x^4
    pop ebp
    ret

; f3= -1/3
f3p:
    push ebp
    mov ebp, esp
    fld1
    fld qword[c3]
    fdivp st1, st0 
    fchs 
    pop ebp
    ret
/*
 * Copyright @ 1984 - 1993   Josef Heinen
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the modified code.  Modifications are to be distributed
 * as patches to released version.
 *
 * This software is provided "as is" without express or implied warranty.
 *
 * Send your comments or suggestions to
 *  J.Heinen@kfa-juelich.de.
 *
 *
 * FACILITY:
 *
 *	GLI Run-Time System
 *
 * ABSTRACT:
 *
 *	This module contains some VMS system services definitions.
 *
 * AUTHOR:
 *
 *	Josef Heinen	10-JAN-1990
 *
 * VERSION:
 *
 *	V1.0
 *
 */


/* System service entry points */

#ifdef VMS
#define SYS_ASSIGN	SYS$ASSIGN
#define SYS_CANCEL	SYS$CANCEL
#define SYS_DASSGN	SYS$DASSGN
#define SYS_QIOW	SYS$QIOW
#endif

#ifdef _WIN32
#define SYS_ASSIGN	sys_assign
#define SYS_CANCEL	sys_cancel
#define SYS_DASSGN	sys_dassgn
#define SYS_QIOW	sys_qiow
#endif

/* System service status codes */

#define SS__NORMAL 1
#define SS__TIMEOUT 556
#define SS__HANGUP 716
#define SS__CANCEL 2096
#define SS__CONTROLC 1617
#define SS__RESIGNAL 2328

/* QIO function codes */

#define IO__WRITEVBLK 48
#define IO__READVBLK 49

/* QIO function modifiers */

#define IO_M_NOECHO 64
#define IO_M_TIMED 128
#define IO_M_NOFORMAT 256
#define IO_M_NOFILTR 512

/* Carriage control information */

#ifdef VMS
#define PREFIX "\r\n"
#define POSTFIX "\r"
#else
#define PREFIX "\r"
#define POSTFIX "\r\n"
#endif

/* RMS status codes */

#define RMS__NORMAL	0x00010001
#define RMS__FNF	0x00018292
#define RMS__ACC	0x0001C002
#define RMS__FND	0x0001C02a

/* Status message definitions */

#define STS_M_SEVERITY	0x00000007
#define STS_M_COND_ID	0x0ffffff8
#define STS_M_CONTROL	0xf0000000
#define STS_M_SUCCESS	0x00000001
#define STS_M_MSG_NO	0x0000fff8
#define STS_M_CODE	0x00007ff8
#define STS_M_FAC_SP	0x00008000
#define STS_M_CUST_DEF	0x08000000
#define STS_M_INHIB_MSG	0x10000000
#define STS_M_FAC_NO	0x0fff0000

#define STS_K_WARNING	0
#define STS_K_SUCCESS	1
#define STS_K_ERROR	2
#define STS_K_INFO	3
#define STS_K_SEVERE	4
#define STS_K_FATAL	4

#define STS_S_CODE      0x0C
#define STS_S_COND_ID   0x19
#define STS_S_CONTROL   0x04
#define STS_S_FAC_NO    0x0C
#define STS_S_MSG_NO    0x0D
#define STS_S_SEVERITY  0x03

#define STS_V_CODE      0x03
#define STS_V_COND_ID   0x03
#define STS_V_CONTROL   0x1C
#define STS_V_CUST_DEF  0x1B
#define STS_V_FAC_NO    0x10
#define STS_V_FAC_SP    0x0F
#define STS_V_INHIB_MSG 0x1C
#define STS_V_MSG_NO    0x03
#define STS_V_SEVERITY  0x00
#define STS_V_SUCCESS   0x00

#define STS_M_MSG_TEXT	(1<<0)
#define STS_M_MSG_ID	(1<<1)
#define STS_M_MSG_SEV	(1<<2)
#define STS_M_MSG_FAC	(1<<3)

/* Define MACROS to extract individual fields from a status value */

#define STATUS_CODE(code)	(((code) & STS_M_CODE) >> STS_V_CODE)
#define STATUS_COND_ID(code)	(((code) & STS_M_COND_ID) >> STS_V_COND_ID)
#define STATUS_CONTROL(code)	(((code) & STS_M_CONTROL) >> STS_V_CONTROL)
#define STATUS_CUST_DEF(code)	(((code) & STS_M_CUST_DEF) >> STS_V_CUST_DEF)
#define STATUS_FAC_NO(code)	(((code) & STS_M_FAC_NO) >> STS_V_FAC_NO)
#define STATUS_FAC_SP(code)	(((code) & STS_M_FAC_SP) >> STS_V_FAC_SP)
#define STATUS_INHIB_MSG(code)	(((code) & STS_M_INHIB_MSG) >> STS_V_INHIB_MSG)
#define STATUS_MSG_NO(code)	(((code) & STS_M_MSG_NO) >> STS_V_MSG_NO)
#define STATUS_SEVERITY(code)	(((code) & STS_M_SEVERITY) >> STS_V_SEVERITY)
#define STATUS_SUCCESS(code)	(((code) & STS_M_SUCCESS) >> STS_V_SUCCESS)

typedef struct {
    int arg_count;
    int name;
    int additional[250];
} CHF_R_SIGNAL_ARGS;

typedef struct {
    int arg_count;
    int *frame;
    int depth;
    int savr0, savr1;
} CHF_R_MECH_ARGS;
 

void sys_authorize(int nusers, char **users);
void sys_listen(unsigned short);
void sys_receive(void);
int sys_socket(void);
void sys_shutdown(void);
int sys_system(char *);
int SYS_ASSIGN(char *, unsigned int *, int, int);
int SYS_DASSGN(int);
int SYS_QIOW(int, int, int, short *, int (*)(), int,
    unsigned char *, int, int, int, char *, int);
int SYS_CANCEL(int);
void establish(int (*)());
void raise_exception(unsigned int, int, ...);
int get_status_text (int, unsigned int, char *);

/*
 * Copyright (C) 2007 Jeremey Murphy & Evan Thomas
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */


/************************************************
   This code implements the NEURON Exp2Syn model
**************************************************/


#include  <parplex/ndl.h>

/*------------------------------
  Exp2Syn - the Neuron synapse
-------------------------------*/
typedef struct {
    DYNAMICS
  double Gmax;    /* The maximum conductance */
  double Er;      /* The driving potential */
  double tau1, tau2;
}  Exp2Syn;


static double Exp2Syn_current(Exp2Syn *self, double t) {
  double v    = GETEM_DYN(self, 0);
  double Gmax = self->Gmax;
  double Er   = self->Er;
  double A    = GETSTATE_DYN(self, 0);
  double B    = GETSTATE_DYN(self, 1);
  return -Gmax * (B-A) * (v-Er);
}

static PyMemberDef Exp2Syn_members[] = {
    {"Gmax", T_DOUBLE, offsetof(Exp2Syn, Gmax), 0, "maximum conductance"},
    {"Er", T_DOUBLE, offsetof(Exp2Syn, Er), 0, "reversal potential"},
    {"tau1", T_DOUBLE, offsetof(Exp2Syn, tau1), 0, "rise time"},
    {"tau2", T_DOUBLE, offsetof(Exp2Syn, tau2), 0, "decay time"},
    {NULL},
};

static void update(Exp2Syn *self, double strength) {
  double A = GETSTATE_DYN(self, 0);
  double B = GETSTATE_DYN(self, 1);
  double factor, tp;
  double tau1 = self->tau1;
  double tau2 = self->tau2;

  if (tau1/tau2 > .9999) { 
    tau1 = .9999*tau2; 
  } 
  tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1);
  factor = -exp(-tp/tau1) + exp(-tp/tau2);
  factor = strength/factor;
 
  A += factor;
  B += factor;
 
  SETSTATE_DYN(self, 0, A);
  SETSTATE_DYN(self, 1, B);
 }

static void Exp2Syn_accepter(Exp2Syn *self, Synapse *s,
			     double strength, int window_id) {
  update(self, strength);
}

static void Exp2Syn_enq(Exp2Syn *self, double starttime,
			double strength) {
  update(self, strength);
}

static void Exp2Syn_derivs(Exp2Syn *self, double t) {
  SETDERIV_DYN(self, 0, -GETSTATE_DYN(self, 0)/self->tau1);
  SETDERIV_DYN(self, 1, -GETSTATE_DYN(self, 1)/self->tau2);
}

static int Exp2SynInit(Exp2Syn *self) {
  SETSTATE_DYN(self, 0, 0);
  SETSTATE_DYN(self, 1, 0);
  return 0;
}

static void HIupdate(Exp2Syn *self, double t, double h) {
  double A = GETSTATE_DYN(self, 0);
  double B = GETSTATE_DYN(self, 1);
  double tau1 = 2*self->tau1;
  double tau2 = 2*self->tau2;

  A *= (tau1-h)/(tau1+h);
  B *= (tau2-h)/(tau2+h);

  SETSTATE_DYN(self, 0, A);
  SETSTATE_DYN(self, 1, B);
}

static void HIcurrent(Exp2Syn *self, double t, double *I, double *Er) {
  double A = GETSTATE_DYN(self, 0);
  double B = GETSTATE_DYN(self, 1);
  *I  = self->Gmax * (B-A);
  *Er = self->Er;
}

static char *Exp2SynStateVars[] = {"A", "B"};
static char *Exp2SynDerivVars[] = {"dAdt", "dBdt"};
static char *Exp2SynTraceVars[] = {"ATrace", "BTrace"};
DynamicsDescriptor Exp2SynDescriptor = {
    "Exp2Syn",
    "Exp2Syn: The neuron synapse model",
    Exp2Syn_members,
    0,
    2,
    Exp2SynStateVars,
    Exp2SynDerivVars,
    Exp2SynTraceVars,
    (derivsfcn*)Exp2Syn_derivs,  /* Derivs */
    (currentfcn*)Exp2Syn_current, /* current */
    (accepterfcn*)Exp2Syn_accepter,  /* accepter */
    (enqfcn*)Exp2Syn_enq,  /* enq */
    0,
    sizeof(Exp2Syn),
    0,
    (userinitfcn*)Exp2SynInit,
    0,
    0,
    0,
    0,
    0,
    0,
    HIdynamics,
    0,
    (HIupdatefcn*)HIupdate,
    (HIcurrentfcn*)HIcurrent    
};
REGISTER_DESCRIPTOR(Exp2SynDescriptor)
  
static DynamicsDescriptor *userDynamics[] = {
  &Exp2SynDescriptor
};

static initproc LuserInitDynamics[] = {
  initExp2SynDescriptor
};

MAKE_P3_MODEL(Exp2Syn, "The NEURON PSC model")

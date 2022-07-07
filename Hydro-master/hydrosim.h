#ifndef _HYDROSIM_H_
#define _HYDROSIM_H_

#include "LinkList.h"



void Simulation(double dt,LinkList<2> &linklist);

void vSimulation(double dt,LinkList<2> &linklist);

void svSimulation(double dt,LinkList<2> &linklist);
void BSQSimulation(double dt,LinkList<2> &linklist);

void Simulation(double dt,LinkList<3> &linklist);


template <int D>
void idealhydro3(LinkList<D> &linklist);

template <int D>
void vischydro3(LinkList<D>  &linklist);

template <int D>
void shear(LinkList<D>  &linklist);
template <int D>
void BSQshear(LinkList<D>  &linklist);


#endif

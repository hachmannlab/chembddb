import numpy as np
import sys
import json

def fetch_unit_list(cur):
    #TODO Fetch from db and parse as a dict of dicts
    cur.execute('USE unit_list_chembddb;')
    cur.execute('SELECT unit_str from MAIN where id=1;')
    unit_list = cur.fetchall()
    print(unit_list)
    unit_list = json.loads(unit_list[0][0])
    return unit_list

def create_unit_list(cur,con):

    property_list = {
    'density':{'g/cm^3 (default)': 1.0, 'kg/m^3': 0.001, 'lb/ft^3': 62.43, 'lb/gal' : 8.35, 'kg/litre': 1.0},
    'energy':{'eV (default)': 1.0, 'Eh': 0.036749, 'J':1.6022e-19, 'Cal/mol':23061.0, 'J/mol':96487.0,'kJ/mol':96.487,'kJ':1.6022e-22,'kCal/mol':23.061, 'Cal':3.8293E-20,'kCal':3.83E-23},
    'polarizability':{'Bohr^3 (default)': 1.0, 'Angstrom': 0.14819},
    'solubility parameters':{'Cal^1/2 cm^-3/2 (default)': 1.0, 'J^1/2 m^-3/2': 2.05 * 10**3, 'MPa^1/2': 2.045},
    'dipole moment':{'Debye (default)': 1.0, 'C-m': 3.34 * 10**30, 'esu-cm': 1 * 10**-18, 'au': 0.3935},
    'quadrupole moment':{'C/m^2 (default)': 1.0, 'C/cm^2': 10000},
    'ratio': {'NA':1.0}
    }
    
    cur.execute('CREATE DATABASE unit_list_chembddb;')
    cur.execute('CREATE TABLE unit_list_chembddb.Main(`id` INT NOT NULL AUTO_INCREMENT, `unit_str` VARCHAR(10000) DEFAULT \'NONE\', PRIMARY KEY (`id`));')
    insert_unit_list(cur,con,property_list)
    return property_list

def insert_unit_list(cur,con,property_list):
    fordb = json.dumps(property_list)
    cur.execute('USE unit_list_chembddb;')
    cur.execute('INSERT INTO Main(id, unit_str) VALUE(1,%s);',[fordb])
    con.commit()

def unit_converter(prop, input_unit, input_val, output_unit):
    f1 = False
    if type(input_value) == float or type(input_value) == list:
        f1 = True
        
    if f1 == False:
        print("Please enter a valid data type of float or list")
        sys.exit()
    
    
    output_value = []
    density = {'g/cm^3': 1.0, 'kg/m^3': 1000.0, 'lb/ft^3': 62.43, 'lb/gal' : 8.35, 'kg/litre': 1.0}
    #dynamic_viscosity = {'kg_per_meter_second': 1.0, 'Poise': 0.1}
    #kinematic_viscosity = {'Stoke': 1.0, 'm^2_per_second': 1*10**-3}
    #mw = {'gram_per_mol': 1.0, 'kg_per_mol': 0.001, 'lb_per_mol' : 0.002205}
    #vdwr = {'pm': 1.0, 'angstrom': 0.01, 'm' : 1* 10**-12, 'cm': 1*10**-10, 'mm': 1*10**-9}
    #ir = {'pm': 1.0, 'angstrom': 0.01, 'm' : 1* 10**-12, 'cm': 1*10**-10, 'mm': 1*10**-9}
    en = {'eV': 1.0, 'Eh': 0.036749, 'J':1.6022e-19, 'Cal/mol':23061.0, 'J/mol':96487.0,'kJ/mol':96.487,'kJ':1.6022e-22,'kCal/mol':23.061, 'Cal':3.8293E-20,'kCal':3.83E-23}
    #pres = {'bar': 1.0, 'Pa': 1* 10**5, 'mmHg' : 750.062, 'psi': 14.5, 'atm': 0.987, 'torr': 750.062}
    #am = {'Da': 1.0, 'kg': 1.66 * 10**-27, 'amu' : 182.888, 'MeV_per_csquare': 931.494}
    pv = {'Bohr^3': 1.0, 'Angstrom': 0.14819}
    sp = {'Cal^1/2 cm^-3/2': 1.0, 'J^1/2 m^-3/2': 2.05 * 10**3, 'MPa^1/2': 2.045}
    dm = {'Debye': 1.0, 'C-m': 3.34 * 10**30, 'esu-cm': 1 * 10**-18, 'au': 0.3935}
    qm = {'C/m^2': 1.0, 'C/cm^2': 10000}
    
    # Unit conversion for density
    if prop.lower() in ['density','number density']:
        #output_value = []
    
        if type(input_value) == list:
            input_val = list(input_value)
            for value in input_val:
                if input_unit == 'gram_per_cm^3':
                    t1 = density[output_unit]
                    cal = t1 * value
                    output_value.append(cal)
                else:
                    #Converting to reference units
                    t1 = (1 / density[input_unit]) * value
                    t2 = t1 * density[output_unit]
                    output_value.append(t2)
                
    
        if type(input_value) == float:
            if input_unit == 'gram_per_cm^3':
                t1 = density[output_unit]
                cal = t1 * input_value
                output_value.append(cal)
            else:
                t1 = (1/ density[input_unit]) * input_value
                t2 = t1 * density[output_unit]
                output_value.append(t2)
            
    # Unit conversion for dynamic viscosity           
    #if prop.lower() == 'dynamic_viscosity':
        
    #    #output_value = []
        
    #    if type(input_value) == list:
    #        input_val = list(input_value)
    #        for value in input_val:
    #            if input_unit == 'kg_per_meter_second':
    #               t1 = dynamic_viscosity[output_unit]
    #               cal = t1 * value
    #               output_value.append(cal)
    #            else:
    #                t1 = (1 / dyanmic_viscosity[input_unit]) * value
    #                t2 = t1 * dynamic_viscosity[output_unit]
    #                output_value.append(t2)
                
    
    #    if type(input_value) == float:
    #        if input_unit == 'kg_per_meter_second':
    #            t1 = dynamic_viscosity[output_unit]
    #            cal = t1 * input_value
    #            output_value.append(cal)
    #        else:
    #            t1 = (1/dynamic_viscosity[input_unit]) * input_value
    #            t2 = t1 * dynamic_viscosity[output_unit]
    #            output_value.append(t2)        
    
    # Unit conversion for kinematic viscosity           
    #if prop == 'kinematic_viscosity':
        
    #    #output_value = []
        
    #    if type(input_value) == list:
    #        input_val = list(input_value)
    #        for value in input_val:
    #            if input_unit == 'Stoke':
    #               t1 = kinematic_viscosity[output_unit]
    #               cal = t1 * value
    #               output_value.append(cal)
    #            else:
    #                t1 = (1 /kinematic_viscosity[input_unit])* value
    #                t2 = t1 * kinematic_viscosity[output_unit]
    #                output_value.append(t2)
                
    
    #    if type(input_value) == float:
    #        if input_unit == 'Stoke':
    #            t1 = Kinematic_viscosity[output_unit]
    #            cal = t1 * input_value
    #            output_value.append(cal)
    #        else:
    #            t1 = (1/kinematic_viscosity[input_unit]) * input_value
    #            t2 = t1 * kinematic_viscosity[output_unit]
    #            output_value.append(t2)
 
    # Unit conversion for Molecular Weight           
    #if prop == 'Molecular Weight':
        
    #    #output_value = []
        
    #    if type(input_value) == list:
    #        input_val = list(input_value)
    #        for value in input_val:
    #            if input_unit == 'gram_per_mol':
    #               t1 = mw[output_unit]
    #               cal = t1 * value
    #               output_value.append(cal)
    #            else:
    #                t1 = (1 / mw[input_unit])* value
    #                t2 = t1 * mw[output_unit]
    #                output_value.append(t2)
                
    
    #    if type(input_value) == float:
    #        if input_unit == 'gram_per_mol':
    #            t1 = mw[output_unit]
    #            cal = t1 * input_value
    #            output_value.append(cal)
    #        else:
    #            t1 = (1/ mw[input_unit])*input_value
    #            t2 = t1 * mw[output_unit]
    #            output_value.append(t2) 
                
    # Unit conversion for Van Der Waals Radius           
    #if prop == 'Van Der Waals Radius':
        
    #    #output_value = []
        
    #    if type(input_value) == list:
    #        input_val = list(input_value)
    #        for value in input_val:
    #            if input_unit == 'pm':
    #               t1 = vdwr[output_unit]
    #               cal = t1 * value
    #               output_value.append(cal)
    #            else:
    #                t1 = (1 /vdwr[input_unit])* value
    #                t2 = t1 * vdwr[output_unit]
    #                output_value.append(t2)
                
    
    #    if type(input_value) == float:
    #        if input_unit == 'pm':
    #            t1 = vdwr[output_unit]
    #            cal = t1 * input_value
    #            output_value.append(cal)
    #        else:
    #            t1 = (1/vdwr[input_unit])* input_value
    #            t2 = t1 * vdwr[output_unit]
    #            output_value.append(t2)
                
    # Unit conversion for Ionic Radius           
    #if prop == 'Ionic Radius':
        
    #    #output_value = []
        
    #    if type(input_value) == list:
    #        input_val = list(input_value)
    #        for value in input_val:
    #            if input_unit == 'pm':
    #               t1 = ir[output_unit]
    #               cal = t1 * value
    #               output_value.append(cal)
    #            else:
    #                t1 = (1/ir[input_unit])*value
    #                t2 = t1 * ir[output_unit]
    #                output_value.append(t2)
                
    
    #    if type(input_value) == float:
    #        if input_unit == 'pm':
    #            t1 = ir[output_unit]
    #            cal = t1 * input_value
    #            output_value.append(cal)
    #        else:
    #            t1 = (1/ir[input_unit])* input_value
    #            t2 = t1 * ir[output_unit]
    #            output_value.append(t2)
                
    # Unit conversion for Energy           
    if prop.lower() == 'energy':
        
        #output_value = []
        
        if type(input_value) == list:
            input_val = list(input_value)
            for value in input_val:
                if input_unit == 'Joules':
                   t1 = en[output_unit]
                   cal = t1 * value
                   output_value.append(cal)
                else:
                    t1 = (1 /en[input_unit]) * value
                    t2 = t1 * en[output_unit]
                    output_value.append(t2)
                
    
        if type(input_value) == float:
            if input_unit == 'Joules':
                t1 = en[output_unit]
                cal = t1 * input_value
                output_value.append(cal)
            else:
                t1 = (1/en[input_unit])* input_value
                t2 = t1 * en[output_unit]
                output_value.append(t2)
                
    # Unit conversion for Pressure           
    if prop == 'Pressure':
        
        #output_value = []
        
        if type(input_value) == list:
            input_val = list(input_value)
            for value in input_val:
                if input_unit == 'bar':
                   t1 = pres[output_unit]
                   cal = t1 * value
                   output_value.append(cal)
                else:
                    t1 = (1 /pres[input_unit])* value
                    t2 = t1 * pres[output_unit]
                    output_value.append(t2)
                
    
        if type(input_value) == float:
            if input_unit == 'bar':
                t1 = pres[output_unit]
                cal = t1 * input_value
                output_value.append(cal)
            else:
                t1 = (1/pres[input_unit]) * input_value
                t2 = t1 * pres[output_unit]
                output_value.append(t2)
                
    # Unit conversion for Atomic Mass           
    if prop == 'Atomic Mass':
        
        #output_value = []
        
        if type(input_value) == list:
            input_val = list(input_value)
            for value in input_val:
                if input_unit == 'Da':
                   t1 = am[output_unit]
                   cal = t1 * value
                   output_value.append(cal)
                else:
                    t1 = (1 /am[input_unit])* value
                    t2 = t1 * am[output_unit]
                    output_value.append(t2)
                
    
        if type(input_value) == float:
            if input_unit == 'Da':
                t1 = am[output_unit]
                cal = t1 * input_value
                output_value.append(cal)
            else:
                t1 = (1/am[input_unit])* input_value
                t2 = t1 * am[output_unit]
                output_value.append(t2)
                
    # Unit conversion for Polarizability Volume           
    if prop.lower() == 'polarizability':
        
        #output_value = []
        
        if type(input_value) == list:
            input_val = list(input_value)
            for value in input_val:
                if input_unit == 'Bohr^3':
                   t1 = pv[output_unit]
                   cal = t1 * value
                   output_value.append(cal)
                else:
                    t1 = (1 /pv[input_unit]) * value
                    t2 = t1 * pv[output_unit]
                    output_value.append(t2)
                
    
        if type(input_value) == float:
            if input_unit == 'Bohr^3':
                t1 = pv[output_unit]
                cal = t1 * input_value
                output_value.append(cal)
            else:
                t1 = (1/pv[input_unit])* input_value
                t2 = t1 * pv[output_unit]
                output_value.append(t2)
                
    # Unit conversion for Solubility Parameter          
    if prop.lower() == 'solubility parameter':
        
        #output_value = []
        
        if type(input_value) == list:
            input_val = list(input_value)
            for value in input_val:
                if input_unit == 'cal^1/2 cm^-3/2':
                   t1 = sp[output_unit]
                   cal = t1 * value
                   output_value.append(cal)
                else:
                    t1 = (1/sp[input_unit]) *value
                    t2 = t1 * sp[output_unit]
                    output_value.append(t2)
                
    
        if type(input_value) == float:
            if input_unit == 'cal^1/2 cm^-3/2':
                t1 = sp[output_unit]
                cal = t1 * input_value
                output_value.append(cal)
            else:
                t1 = (1/sp[input_unit]) *input_value
                t2 = t1 * sp[output_unit]
                output_value.append(t2)


    # Unit conversion for Dipole Moment           
    if prop.lower() == 'dipole moment':
        
        #output_value = []
        
        if type(input_value) == list:
            input_val = list(input_value)
            for value in input_val:
                if input_unit == 'debye':
                   t1 = dm[output_unit]
                   cal = t1 * value
                   output_value.append(cal)
                else:
                    t1 = (1 /dm[input_unit]) * value
                    t2 = t1 * dm[output_unit]
                    output_value.append(t2)
                
    
        if type(input_value) == float:
            if input_unit == 'debye':
                t1 = dm[output_unit]
                cal = t1 * input_value
                output_value.append(cal)
            else:
                t1 = (1/dm[input_unit])* input_value
                t2 = t1 * dm[output_unit]
                output_value.append(t2)


    # Unit conversion for Quadrupole Moment           
    if prop.lower() == 'quadrupole moment':
        
        #output_value = []
        
        if type(input_value) == list:
            input_val = list(input_value)
            for value in input_val:
                if input_unit == 'coulomb_per_m^2':
                   t1 = qm[output_unit]
                   cal = t1 * value
                   output_value.append(cal)
                else:
                    t1 = (1 /qm[input_unit]) * value
                    t2 = t1 * qm[output_unit]
                    output_value.append(t2)
                
    
        if type(input_value) == float:
            if input_unit == 'coulomb_per_m^2':
                t1 = qm[output_unit]
                cal = t1 * input_value
                output_value.append(cal)
            else:
                t1 = (1/qm[input_unit])* input_value
                t2 = t1 * qm[output_unit]
                output_value.append(t2)

                
    # Unit conversion for Temperature          
    if prop.lower() == 'temperature':
        #output_value = []
        
        if type(input_value) == list:
            input_val = list(input_value)
            for value in input_val:
                if input_unit == 'K':
                    if output_unit == 'C':
                        cal = value -273.15
                        output_value.append(cal)
                    if output_unit == 'F':
                        cal = (value - 273.15) * 9/5 + 32
                        output_value.append(cal)
                elif input_unit == 'C':
                    if output_unit == 'K':
                        cal = value + 273.15
                        output_value.append(cal)
                    if output_unit == 'F':
                        cal = (value * 9/5) + 32
                        output_value.append(cal)
                elif input_unit == 'F':
                    if output_unit == 'C':
                        cal = (value - 32) * 5/9
                        output_value.append(cal)
                    if output_unit == 'K':
                        cal = (value - 32) * 5/9 + 273.15
                        output_value.append(cal)
    
        if type(input_value) == float:
            if input_unit == 'K':
                if output_unit == 'C':
                    cal = input_value -273.15
                    output_value.append(cal)
                if output_unit == 'F':
                    cal = (input_value - 273.15) * 9/5 + 32
                    output_value.append(cal)
            elif input_unit == 'C':
                if output_unit == 'K':
                    cal = inpput_value + 273.15
                    output_value.append(cal)
                if output_unit == 'F':
                    cal = (input_value * 9/5) + 32
                    output_value.append(cal)
            elif input_unit == 'F':
                if output_unit == 'C':
                    cal = (input_value - 32) * 5/9
                    output_value.append(cal)
                if output_unit == 'K':
                    cal = (input_value - 32) * 5/9 + 273.15
                    output_value.append(cal)
            

    return output_value
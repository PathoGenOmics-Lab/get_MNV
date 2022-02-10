#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""                                                                                                                                                                                                                                                                                                                                                                                       
                                                                                           .            .       :            ...   ::            ::.  
                                                                                           =:          :=       *+:          **-   %@*         .%@#   
                                                                                         .===.    .   -==.    .:***=.    .  :**:   -@@-.   .. :@@#    
                                                                                         -===-      .====.     =*++*+:      =**     #@#   .  -@@+.    
                                                                                    .   .==:==:    .==-==:  . .**=.:**+.   .**=    .-@@-  ..=@@=..    
                                                                                      ..==-.-== ..:==-.-=-   ::**....+**-. .**:.    .*@@. .*@@:  ..   
                                                                   .****:              -==   ==- :==:  :==    -**     :**=.-**    .  .@@* #@%:        
                                                                 .**.  =%             :==.   .==-==.   .==    +*=   .  .=****=        =@@%@#. .       
                                                                +#:  .**.     #+    .:==: . . -===. ...:==. :.**:..   ...:***-.    ... %@@#.   ..   :@
                                                         .-----%=   +#:       #*     -==       ==       ==:  :**.          =*          :@+          :@
                            :.                           .:::-:    =-         #*               .              ..            .           :           :@
               :=******: -#*-*%           -+*****+.        .%=  .******:      ##====================================================================*@
            .+#+:  :=: :#+   *+        -**-. -+: =%       -%:  =%:            ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
          :**:  =*+-..+#.   :%       =#=  =*+-   =#     .#*  :#=                                     .                                                
        .#*.  =%-   .%=     %-     =%- .*#-     =%.    =%: .*#.                            .::::   -.                                                 
       +#.  .%=    -%: -%  +#    .#+  *#:     -#+    :%+  -%-                            :-:.:::-.                                                    
      +-   =%:    **   %= :@.   =%: :%-:..:=**=     +*.  *#         -%:             .:..:-                                                            
          -%    :%-   =#  +- .=#-  :@. :==-:     .+#:   +#        :#*            ::-:::-:                                                             
      .   %-   +#. : :@.  :+*+:    +#         .=#*:     %=      :#*:               :             -::                                                  
     -@   +*-+#- :#=:@: .%=.    =+ .#*=--=++**=:     *: :%+--=**=.        :::::::::::::::::::::. : .                                                  
     .@:   .:. :** :%:  %=      .%=   .::.  .=**:    *#   .::..-*#:      -::::::::::::::::::: .=                                                      
      :**=--=**=. -%.  #+         -***+++***+:        =#*+=++*+-                                                                                      
         .::.:   +#.  #=                                  .                                                                                           
          -*#= .#+  .%=                                                                                                                               
       -*#=   =%:  :%:                                                                                                                                
    .*#=    :#=   +#.                                                                                                                                 
  .**.    :#+   :%=                                                                                                                                   
  @:   .=#+   :#*                                                                                                                                     
  *+++*=:   -*+.                                                                                                                                      
@-      .-*#=                                                                                                                                         
:**++***+-                                                                                                                                            

Script for annotate vcf file with multiple-nucleotide variants

**IMPORTANT these script works with SNV against a reference, 
insertions and deletions modifiying reading frame are not currently supported**

Command:
    python3 get_mnv.py -v G35894.var.snp.vcf -f MTB_ancestor.fas -g anot_genes.txt
Outputs:
    vcf_name + .MNV.vcf  -> annotated vcfs
    vcf_name + .MNV.tsv  -> tsv file with MNV variants
 
Requires biopython

Contact -> @paururo
"""

__author__ = 'Paula Ruiz Rodriguez'
__credits__ = ['Paula Ruiz Rodriguez', 'PathoGenOmics']
__license__ = 'GPL'
__version__ = '0.1.0'
__maintainer__ = 'Paula Ruiz Rodriguez: @paururo'
__email__ = 'paula.ruiz-rodriguez@uv.es'
__status__ = 'developing'
__lab__ = 'PathoGenOmics, I2SysBio'
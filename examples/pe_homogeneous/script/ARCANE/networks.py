"""Module containing all classes and functions related to kinetic networks"""

# Import statements
import cantera as ct
import numpy as np


class Network:

    def __init__(self, mechanism):
        
        """
        Organize properly info on kinetic networks.


        Created: 11/01/17 [QC]
        Last modified: 01/15/18 [QC]
        """

        # Keep track of which mechanism was used to create network
        self.ctmech = mechanism.ctmech

        # Set up stoichiometric coefficients
        self.nur = ct.Kinetics.reactant_stoich_coeffs(self.ctmech)
        self.nup = ct.Kinetics.product_stoich_coeffs(self.ctmech)

        # Coupling between species
        self.create_deltaSS()

        # Coupling between species and reaction
        # deltaSR[Si, Rj] = 1 if Si participates in reaction Rj
        self.create_deltaSR()

        # Compact lists of indices for species/reactions
        self.inds, self.indr = np.nonzero(self.deltaSR)

        # Compact lists of indices for species/species
        self.indsi, self.indsj = np.nonzero(self.deltaSS)

    def create_deltaSS(self):
        """
        Species/species connectivity matrix
        Matrix of size ns x ns contains 1 values if species i and species j are present in a common reaction

        Parameters
        ----------
        ctmech :
            Cantera Solution object

        Returns
        -------
        type :
            square matrix of dimension number of species with zeros and ones

        """

        # Identify directly linked species in phase
        # Initialize empty list
        self.deltaSS = np.zeros([self.ctmech.n_species, self.ctmech.n_species], 'd')

        # Loop through reactions
        for reac in self.ctmech.reactions():
            # Concatenate reactants and products
            slist = []
            for spec, value in reac.reactants.items():
                slist.append(spec)
            for spec, value in reac.products.items():
                slist.append(spec)

            # For each species in reaction, append the reac name to the right [Si,Sj] entry
            for spec in slist:
                for spec2 in slist:
                    self.deltaSS[self.ctmech.species_index(spec)][self.ctmech.species_index(spec2)] = 1

    def create_deltaSR(self):
        """
        Species/reactions connectivity matrix.      
        
        
        Parameters
        ----------
        Matrix of size ns x nr, contains 1 values if Species i is present in reaction j

        Returns
        -------
        type :
            matrix of dimension (number of species * number of reactions) with zeros and ones

        """

        nu = np.array(abs(self.nur + self.nup))

        nu[nu >= 1] = 1

        self.deltaSR = nu

##### FROM SPECTRUM MODEL CLASS ################################################

    def PlotSpectrum(self):
        # Generate Colours
        Col = GenerateColours()
        plt.figure(figsize=(8,8))
        plt.title('Spectrum')
        plt.xlabel('Wavelength / Å')
        plt.ylabel('Flux')
        plt.plot(self.Wavelengths, self.Fluxes, color=Col[0], lw=.2, label='Spectrum Data')
        plt.xlim(self.Wavelengths[0], self.Wavelengths[-1])
        if self.Errors != None:
            plt.plot(self.Wavelengths, self.Errors, color=Col[1], lw=.2, label='Error')
            plt.legend()
        if self.Metallicity is not None:
            FigureName = f'SpectrumData/Plots/Spectrum[{self.Metallicity}].png'
        else:
            FigureName = 'SpectrumData/Plots/Spectrum.png'
        plt.savefig(FigureName)
        plt.show()

    def WriteOutData(self):
        # Assumes errors within the data
        ColumnNames = ['Wl', 'Flux', 'Errors']
        Data = [self.Wavelengths, self.Fluxes, self.Errors]
        if self.Metallicity is not None:
            Name = f'SpectrumData/Data[{self.Metallicity}].dat'
        else:
            Name = 'SpectrumData/Data.dat'
        OutputData = Table(Data, names=ColumnNames)
        OutputData.write(Name, format='ascii', overwrite=True)

##### FROM SPLINE MODEL CLASS ##################################################

    def PlotSpectrum(self):
        # Generate Colours
        Col = GenerateColours()
        plt.figure(figsize=(8,8))
        plt.title('Spectrum')
        plt.xlabel('Wavelength / Å')
        plt.ylabel('Flux')
        plt.plot(self.Wavelengths, self.Fluxes, color=Col[0], lw=.2, label='Spectrum Data')
        plt.xlim(self.Wavelengths[0], self.Wavelengths[-1])
        if self.Errors != None:
            plt.plot(self.Wavelengths, self.Errors, color=Col[1], lw=.2, label='Error')
            plt.legend()
        if self.Metallicity is not None:
            FigureName = f'SpectrumData/Plots/Spectrum[{self.Metallicity}].png'
        else:
            FigureName = 'SpectrumData/Plots/Spectrum.png'
        plt.savefig(FigureName)
        plt.show()

    def WriteOutData(self):
        # Assumes errors within the data
        ColumnNames = ['Wl', 'Flux', 'Errors']
        Data = [self.Wavelengths, self.Fluxes, self.Errors]
        if self.Metallicity is not None:
            Name = f'SpectrumData/Data[{self.Metallicity}].dat'
        else:
            Name = 'SpectrumData/Data.dat'
        OutputData = Table(Data, names=ColumnNames)
        OutputData.write(Name, format='ascii', overwrite=True)

     def PlotSplineFluxes(self, Windows=None):
         # Select Colours
         Col = GenerateColours()
         plt.figure(figsize=(16,8))
         if self.Metallicity is not None:
             plt.title('S99 Model: Spectrum with Spline [%.2f $Z\\textsubscript{\(\odot\)}$]' % self.Metallicity)
         else:
             plt.title('S99 Model: Spectrum with Spline')
         plt.xlabel('Wavelength / Å')
         plt.ylabel('Flux')

         # If Windows Given, shade the relevant regions
         if Windows is not None:
             WindowData = Table.read(Windows, format='ascii')
             WindowsStart = WindowData['wlmin']
             WindowsStop = WindowData['wlmax']
             for i in range(WindowsStart.shape[0]):
                 plt.axvspan(WindowsStart[i], WindowsStop[i], color='#C0C0C0', alpha=0.2)

         plt.plot(self.Wavelengths, self.Fluxes, color=Col[0], lw=.2, label='Spectrum Data')
         plt.plot(self.Wavelengths, self.SplineFluxes, color=Col[1], lw=.2, label='Spline')
         #plt.plot(self.Wavelengths, self.Errors, color=Col[2], lw=.2, label='Error')
         plt.xlim(self.Wavelengths[0], self.Wavelengths[-1])
         plt.legend()
         if self.Metallicity is not None:
             FigureName = f'Continuum Fitting/Plots/SplinePlot[{self.Metallicity}].png'
         else:
             FigureName = 'TESTS/ModPlot2.png'
         #plt.xlim(WindowsStart[1]-50, WindowsStop[3]+50)
         plt.savefig(FigureName)
         plt.show()

     def PlotRelativeFluxes(self, Windows=None):
         plt.figure(figsize=(8,8))
         if self.Metallicity is not None:
             plt.title('S99 Model: Relative Flux [%.2f $Z\\textsubscript{\(\odot\)}$]' % self.Metallicity)
         else:
             plt.title('S99 Model: Relative Flux')
         plt.xlabel('Wavelength / Å')
         plt.ylabel('Relative Flux')

         # If Windows Given, shade the relevant regions
         if Windows != None:
             WindowData = Table.read(Windows, format='ascii')
             WindowsStart = WindowData['wlmin']
             WindowsStop = WindowData['wlmax']
             for i in range(WindowsStart.shape[0]):
                 plt.axvspan(WindowsStart[i], WindowsStop[i], color='skyblue', alpha=0.2)

         plt.plot(self.Wavelengths, self.RelativeFluxes, color='black', lw=.2, label='RelativeFlux')
         plt.xlim(self.Wavelengths[0], self.Wavelengths[-1])
         plt.hlines(1.0, self.Wavelengths[0], self.Wavelengths[-1], color='black', lw=.1)
         plt.legend()
         if self.Metallicity is not None:
             FigureName = f'Continuum Fitting/Plots/SplineRelFluxPlot[{self.Metallicity.round(2)}].png'
         else:
             FigureName = 'Continuum Fitting/Plots/SplineRelFluxPlot.png'
         plt.ylim(0.5, 1.5)
         plt.savefig(FigureName)
         plt.show()

     def WriteOutSplineData(self):
         ColumnNames = ['Wl', 'Flux', 'SplineFlux', 'RelFlux']
         Data = [self.Wavelengths, self.Fluxes, self.SplineFluxes, self.RelevantFluxes]
         if self.Metallicity is not None:
             Name = f'Continuum Fitting/Data/SplineData[{self.Metallicity.round(2)}].dat'
         else:
             Name = 'Continuum Fitting/Data/SplineData.dat'
         OutputData = Table(Data, names=ColumnNames)
         OutputData.write(Name, format='ascii', overwrite=True)

    def PlotInputVsOutput(InputZ, OutputZ, OutputErrors, WriteOut=False):

        if WriteOut == True:
            # Also will write data to a table
            OutputTable = Table([InputZ, OutputZ, OutputErrors], names=['Z_In, Z_Out, Z_Err'])
            OutputTalbe.write('Chi Square Fitting/Data/InputVsOutputZ.dat', format='ascii', overwrite=True)

        # Generate Plot
        plt.figure(figsize=(8,8))
        plt.title('Comparison of Input vs Output Metallicites')
        plt.scatter(InputZ, OutputZ, c='black', s=0.4)
        plt.errorbar(InputZ, OutputZ, yerr=OutputErrors, fmt='none', elinewidth=.5, ecolor='black', capsize=4)
        plt.errorbar(InputZ, OutputZ, yerr=(2*OutputErrors), fmt='none', elinewidth=.25, ecolor='red', capsize=2)
        plt.plot(InputZ, OutputZ, linestyle='--', color='black', lw=.2)
        plt.xlabel("Input Metallicity / $Z\\textsubscript{\(\odot\)}$")
        plt.ylabel('Output Metallicity / $Z\\textsubscript{\(\odot\)}$')
        plt.xlim(0.0, 3.0)
        plt.ylim(0.0, 3.0)
        plt.savefig('Chi Square Fitting/Plots/InputVsOutputZ.png')
        plt.show()

        # Plot Difference Plot
        DeltaZ = OutputZ - InputZ
        plt.figure(figsize=(8,8))
        plt.title('Comparison of Input vs Output Metallicites')
        plt.scatter(InputZ, DeltaZ, c='black', s=0.4)
        plt.errorbar(InputZ, DeltaZ, yerr=OutputErrors, fmt='none', elinewidth=.5, ecolor='black', capsize=4)
        plt.errorbar(InputZ, DeltaZ, yerr=(2*OutputErrors), fmt='none', elinewidth=.25, ecolor='red', capsize=2)
        plt.hlines(0.0, 0.0, 3.0, color='black', lw=.1)
        plt.xlabel("Input Metallicity / $Z\\textsubscript{\(\odot\)}$")
        plt.ylabel('Difference in Metallicity / $Z\\textsubscript{\(\odot\)}$')
        plt.xlim(0.0, 3.0)
        plt.savefig('Chi Square Fitting/Plots/InputVsOutputZ_Comparison.png')
        plt.show()

    def PlotComparison(ModelArray, Windows=None):
        # Set Figure Parameters
        plt.figure(figsize=(16,8))
        plt.title('S99 Model Comparison: Relative Flux')
        plt.xlabel('Wavelength / Å')
        plt.ylabel('Relative Flux  + Constant')
        OffSet = 0.0

        # If Windows Given, shade the relevant regions
        if Windows != None:
            WindowData = Table.read(Windows, format='ascii')
            WindowsStart = WindowData['wlmin']
            WindowsStop = WindowData['wlmax']
            for i in range(WindowsStart.shape[0]):
                plt.axvspan(WindowsStart[i], WindowsStop[i], color='skyblue', alpha=0.2)

        # Plot Relevant Data
        for Model in ModelArray:
            ModelFluxes = Model.RelativeFluxes + np.repeat(OffSet, Model.RelativeFluxes.shape[0])
            plt.plot(Model.Wavelengths, ModelFluxes, color='black', lw=.2)
            if Model.Metallicity is not None:
                Label = "%.2f $Z\\textsubscript{\(\odot\)}$" % Model.Metallicity
                plt.text(Model.Wavelengths[-1]+10, OffSet+0.98, Label)
            plt.hlines(OffSet + 1.0, Model.Wavelengths[0], Model.Wavelengths[-1], color='black', lw=.1)
            OffSet += 1.0

        # Set Figure Parameters
        plt.xlim(ModelArray[0].Wavelengths[0], ModelArray[0].Wavelengths[-1])
        plt.ylim(0.0, OffSet+1.0)
        plt.savefig('Continuum Fitting/Plots/ModelComparisonPlot.png')
        plt.show()

    def ChiSquare(self, ModelFluxes, WindowsPath):

        # Determine the Chi Square Windows
        ContinuumWindows = Table.read(WindowsPath, format='ascii')
        WindowWlMins = ContinuumWindows['wlmin']
        WindowWlMaxs = ContinuumWindows['wlmax']

        # Set Loop parameters and arrays
        WindowNumber = WindowWlMins.shape[0]
        ChiSquareArray = np.empty(0)
        Count = 0

        # Iterate over all pixels in windows
        for window in range(WindowNumber):
            for i in range(self.Wavelengths.shape[0]):
                wl = self.Wavelengths[i]
                if wl < WindowWlMins[window] or wl > WindowWlMaxs[window]:
                    continue
                else:
                    Count += 1
                    ModelValue = ModelFluxes[i]
                    DataValue = self.RelativeFluxes[i]
                    ErrorValue = self.Errors[i]
                    ChiSquare = ((ModelValue-DataValue)**2)/(ErrorValue**2)
                    ChiSquareArray = np.append(ChiSquareArray, ChiSquare)
        # Sum the Chi array
        ChiValue = np.sum(ChiSquareArray)
        # Caluclate the Reduced Chi
        ReducedChi = ChiValue / Count
        # Return Value to the Main
        return ReducedChi

    def GenerateSteidelMask(Wavelengths):
        ContinuumWindows = Table.read('Parameter Files/Windows/SteidelWindows.txt', format='ascii')
        WindowWlMins = ContinuumWindows['wlmin']
        WindowWlMaxs = ContinuumWindows['wlmax']

        Mask = np.array([np.any((WindowWlMins <= wl) & (wl <= WindowWlMaxs)) for wl in Wavelengths])

        return Mask

    def MonteCarloError(StackData, Sweeps, ModelArray, WindowsPath):

        ##### NEEDS FIXED
        # Determine Chi Square without perturbation
        Metallicities = np.empty(Sweeps)

        # Iterate over the sweeps
        for i in range(Sweeps):
            # Chi Array for individual sweeps
            ChiArray = np.empty(ModelArray.shape[0])
            # Perturb Fluxes if not initial sweep
            OriginalFluxes = StackData.RelativeFluxes
            if i == 0:
                StackData.RelativeFluxes = StackData.RelativeFluxes
            else:
                StackData.RelativeFluxes = SpectrumModel.PerturbFluxes(StackData.RelativeFluxes, StackData.Errors)
            # Iterate over all models getting the Chi Square Values
            for j in range(ModelArray.shape[0]):
                ChiArray[i] = StackData.ChiSquare(ModelArray[j].RelativeFluxes, WindowsPath)
            # Determinine Metallicity at the Minimum Chi
            Metallicities = np.append(Metallicities, 10**ModelArray[np.argmin(ChiArray)].Metallicity)
            print("[ Completion : " + str(round(100*(i)/(Sweeps),1)) + '%]', end= '\r')

        # Return Output Metallicity array to main
        return Metallicities

    def GetEstimates(Metallicities):
        # Determine Values from the Metallicities Array
        Mean = np.mean(Metallicities)
        Error = np.std(Metallicities)
        return Mean, Error

    def PlotMCDistribution(Metallicities, Bins, Spacing=0.1):
        plt.figure(figsize=(8,8))
        plt.hist(Metallicities, Bins, density=True, color='black')
        plt.xlabel('Metallicity')
        plt.ylabel('Normalised Frequency')
        plt.title('Distribution of Metallicites')
        plt.savefig('Chi Squared Fitting/Plots/DistributionOfZ.png')
        plt.show()
    def PlotBestFitModel(Data, Z, ZErr, ModelArray, Windows=False, TestIndex=0, TrueZ = None):

        # Set colours and opacity
        CArr = GenerateColours()
        A = 0.5

        # Define Model Values on the offchance a model has an equivalent metallicity
        ZArray = np.array([0.001, 0.002, 0.008, 0.014, 0.040])
        ZArraySolar = (ZArray/0.0142).round(2)

        ZArr = np.array([Z - ZErr[0], Z, Z + ZErr[1]])
        ZMods = np.empty(shape=ZArr.shape[0], dtype=object)

        # Interpolate Models to correct Metallicity
        for i in range(ZArr.shape[0]):
            Z = ZArr[i]
            if Z in ZArraySolar:
                Model = ModelArray[np.where(ZArraySolar==Z)[0][0]]
            else:
                Model = SpectrumModel.InterpolateModels(Z, ModelArray, Linear=True)
            # Convolve to correct Resolution
            Model.ConvolveModel(FWHM=3.0, PixelRes=0.4)
            # Determine Relevant Continuum Spline
            SplineModel = Model.DetermineContinuumSpline('Parameter Files/Windows/RixWindows.txt')
            # Resample onto Wavelength Grid
            SplineModel.ResampleRelativeFluxToGrid(Data.Wavelengths)
            ZMods[i] = SplineModel

        fig = plt.figure(figsize=(16,8))
        spec = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[4,1])

        ax0 = fig.add_subplot(spec[0])
        ax1 = fig.add_subplot(spec[1])

        if TrueZ is not None:
            ZLabel = "Data: %.3f $Z\\textsubscript{\(\odot\)}$" % TrueZ
        else:
            ZLabel = "Data"

        # Plot Data
        ax0.plot(Data.Wavelengths, Data.RelativeFluxes, color=CArr[1], lw=.5, label=ZLabel)
        # Plot best fitting model
        ax0.plot(ZMods[1].Wavelengths, ZMods[1].RelativeFluxes, color=CArr[2], lw=.5, label='Best Fit Model: %.3f $Z\\textsubscript{\(\odot\)}$' % ZMods[1].Metallicity)
        # Plot Errors
        #MinDataFlux = Data.RelativeFluxes - Data.Errors
        #MaxDataFlux = Data.RelativeFluxes + Data.Errors
        #ax0.fill_between(Data.Wavelengths, MinDataFlux, MaxDataFlux, alpha=A, color=CArr[1], label='Data Error')
        ax1.plot(Data.Wavelengths, Data.Errors, color=CArr[4], lw=.5, label='Spectral Error')
        #plt.fill_between(Data.Wavelengths, ZMods[2].RelativeFluxes, ZMods[0].RelativeFluxes, alpha=0.7, color='#C0C0C0')
        MinFlux = np.minimum(ZMods[2].RelativeFluxes, ZMods[0].RelativeFluxes)
        MaxFlux = np.maximum(ZMods[2].RelativeFluxes, ZMods[0].RelativeFluxes)
        ax0.fill_between(Data.Wavelengths, MinFlux, MaxFlux, alpha=A, color=CArr[2], label='Best Fit Model Error')
        # Plot hline
        ax0.hlines(1.0, Data.Wavelengths[0], Data.Wavelengths[-1], color='black', lw=.1)
        # Plot Parameters
        ax1.set_xlabel('Wavelength / Å')
        ax0.set_ylabel('Continuum Normalised Flux')
        ax1.set_ylabel('Continuum Normalised Error')
        ax0.set_title(f'Plot of Best Fitting Model for Test Spectrum {TestIndex}')
        ax0.legend(loc='lower right')
        ax1.legend(loc='lower right')
        ax0.set_xlim(Data.Wavelengths[0], Data.Wavelengths[-1])
        ax1.set_xlim(Data.Wavelengths[0], Data.Wavelengths[-1])
        plt.subplots_adjust(wspace=0, hspace=0)
        #ax0.set_ylim(0.6, 1.1)

        # If Windows Given, shade the relevant regions
        if Windows is not None:
            WindowData = Table.read(Windows, format='ascii')
            WindowsStart = WindowData['wlmin']
            WindowsStop = WindowData['wlmax']
            for i in range(WindowsStart.shape[0]-1):
                if i == 0:
                    ax0.axvspan(Data.Wavelengths[0], WindowsStart[i], color='#C0C0C0', alpha=0.2)
                    ax1.axvspan(Data.Wavelengths[0], WindowsStart[i], color='#C0C0C0', alpha=0.2)
                elif i == WindowsStart.shape[0]-1:
                    ax0.axvspan(WindowsStop[i], Data.Wavelengths[-1], color='#C0C0C0', alpha=0.2)
                    ax1.axvspan(WindowsStop[i], Data.Wavelengths[-1], color='#C0C0C0', alpha=0.2)
                else:
                    ax0.axvspan(WindowsStop[i], WindowsStart[i+1], color='#C0C0C0', alpha=0.2)
                    ax1.axvspan(WindowsStop[i], WindowsStart[i+1], color='#C0C0C0', alpha=0.2)

        # Save and display Figure
        plt.savefig(f'Fergus_Data/Plots/BFM/m{TestIndex}.png')

##### END SPLINE MODEL #########################################################                         

import numpy as np
from astropy import units
import healpy as hp
from scipy.interpolate import interp1d

from .component_data import HASLAM_FILEPATH
from .base_skymodel import BaseSkyModel
from .base_observer import BaseObserver


class HaslamSkyModel(BaseSkyModel):
    """ Haslam destriped, desourced sky model """

    def __init__(self,  freq_unit='MHz', spectral_index=-2.6):
        """ Global sky model (GSM) class for generating sky models.

        Parameters
        ----------
        freq_unit (str): Frequency unit to use, defaults to MHz
        spectral_index (float): Spectral index to use for calculations.

        Notes
        -----
        This is a crude model; the sky's spectral index changes with sidereal time
        and a singular spectral index is not realistic. Here are some measured values
        at different frequencies:

            Frequency        Spectral Index        Reference     Notes
            50--100 MHz     si = -2.56 +/- 0.03    Mozdzen+19    LST 0--12 hr
                                  2.46 +/- 0.01                  LST 18.2 hr
            90--190 MHz     si = -2.61 +/- 0.01    Mozdzen+16    LST 0--12 hr
                                 −2.50 +/- 0.02                  LST 17.7 hr
            0.4--4.7 GHz 	si = −2.85 +/- 0.05    Dickinson+19
            0.4--22.8 GHz 	si = −2.88 +/- 0.03    Dickinson+19
            4.7--22.8 GHz 	si = −2.91 +/- 0.04    Dickinson+19
            22.8--44 GHz 	si = −2.85 +/- 0.14    Dickinson+19

        References
        ----------
        Mozdzen et al (2016), EDGES high-band, https://doi.org/10.1093/mnras/stw2696
        Mozdzen et al (2019), EDGES low-band, https://doi.org/10.1093/mnras/sty3410
        Dickinson et al (2019), C-BASS experiment, https://doi.org/10.1093/mnras/stz522

        """
        data_unit = 'K'
        basemap = 'Haslam'
        super(HaslamSkyModel, self).__init__('Haslam', HASLAM_FILEPATH, freq_unit, data_unit, basemap)
        self.spectral_index = spectral_index
        self.data = hp.read_map(self.fits, verbose=False, dtype=np.float64)

    def generate(self, freqs):
        """ Generate a global sky model at a given frequency or frequencies

        Parameters
        ----------
        freqs: float or np.array
            Frequency for which to return GSM model

        Returns
        -------
        gsm: np.array
            Global sky model in healpix format, with NSIDE=256. Output map
            is in galactic coordinates, and in antenna temperature units (K).

        """
        # convert frequency values into Hz
        freqs = np.array(freqs) * units.Unit(self.freq_unit)
        freqs_mhz = freqs.to('MHz').value

        if isinstance(freqs_mhz, float):
            freqs_mhz = np.array([freqs_mhz])

        map_out = self.data * (freqs_mhz / 408.0) ** (self.spectral_index)

        self.generated_map_data = map_out
        self.generated_map_freqs = freqs
        return map_out


class HaslamObserver(BaseObserver):
    def __init__(self):
        """ Initialize the Observer object.

        Calls ephem.Observer.__init__ function and adds on gsm
        """
        super(HaslamObserver, self).__init__(gsm=HaslamSkyModel)
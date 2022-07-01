class TimeSlicing:
    def which(self, cosmo, tau):
        raise NotImplementedError


class TimeSlicingRegion(TimeSlicing):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def which(self, cosmo, tau):
        return (tau >= self.start) & (tau < self.end)


class TimeSlicingReco(TimeSlicing):
    def __init__(self, end):
        """
        end: fraction of tau_rec at which to end
        """
        self.end = end

    def which(self, cosmo, tau):
        tau_rec = cosmo.get_current_derived_parameters(["tau_rec"])["tau_rec"]
        return tau < self.end * tau_rec


class TimeSlicingReio(TimeSlicing):
    def __init__(self, start):
        """
        start: fraction of tau_reio at which to start
        """
        self.start = start

    def which(self, cosmo, tau):
        tau_rec = cosmo.get_current_derived_parameters(["tau_rec"])["tau_rec"]
        z_reio = cosmo.get_current_derived_parameters(['z_reio'])['z_reio']
        tau_reio = cosmo.tau_of_z(z_reio)
        return tau >= self.start * tau_reio


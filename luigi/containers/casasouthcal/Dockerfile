FROM kernsuite/casa:4.7
RUN rm -rf \
    /casa-release-4.7.0-el6/data/nrao/VLA/standards/fluxcalibrator.data \
    /casa-release-4.7.0-el6/data/nrao/VLA/standards/PerleyButler2013Coeffs
COPY fluxcalibrator.data /casa-release-4.7.0-el6/data/nrao/VLA/standards/fluxcalibrator.data
COPY PerleyButler2013Coeffs /casa-release-4.7.0-el6/data/nrao/VLA/standards/PerleyButler2013Coeffs


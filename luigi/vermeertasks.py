from kliko.luigi_util import KlikoTask


class DownloadTask(KlikoTask):
    @classmethod
    def image_name(cls):
        return "vermeerkat/downobs:0.1"


class H5tomsTask(KlikoTask):
    @classmethod
    def image_name(cls):
        return "vermeerkat/h5toms:0.1"

    def requires(self):
        return DownloadTask()


class RfiMaskerTask(KlikoTask):
    @classmethod
    def image_name(cls):
        return "vermeerkat/rfimasker:0.1"

    def requires(self):
        return H5tomsTask()


class AutoFlaggerTask(KlikoTask):
    @classmethod
    def image_name(cls):
        return "vermeerkat/autoflagger:0.1"

    def requires(self):
        return RfiMaskerTask()


class WscleanTask(KlikoTask):
    @classmethod
    def image_name(cls):
        return "vermeerkat/wsclean:0.1"

    def requires(self):
        return AutoFlaggerTask()

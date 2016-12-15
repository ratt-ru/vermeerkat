from kliko.luigi_util import KlikoTask


class SimmsTask(KlikoTask):
    @classmethod
    def imagename(cls):
        return "vermeerkat/h5toms:0.1"


class RfiMaskerTask(KlikoTask):
    @classmethod
    def imagename(cls):
        return "vermeerkat/rfimasker:0.1"

    def requires(self):
        return SimmsTask()


class AutoFlaggerTask(KlikoTask):
    @classmethod
    def imagename(cls):
        return "vermeerkat/autoflagger:0.1"

    def requires(self):
        return RfiMaskerTask()


class WscleanTask(KlikoTask):
    @classmethod
    def imagename(cls):
        return "vermeerkat/wsclean:0.1"

    def requires(self):
        return AutoFlaggerTask()

from .adapters import SourceFileAdapter, CLASSAdapter
from .processors import InputPostProcessor, TargetPostProcessor


class InputProvider:
    def get_inputs(self, tau=None):
        pass


class OutputProvider:
    def get_outputs(self, tau=None):
        pass


class AdapterProvider(InputProvider, OutputProvider):

    def __init__(self, adapter):
        self.adapter = adapter

        # After the inputs have been obtained
        # the first time, they will be cached
        # in this field
        self._adapter_inputs = None

        self.processor = InputPostProcessor()

    def get_adapter_inputs(self):
        if self._adapter_inputs is None:
            self._adapter_inputs = self.adapter.get_inputs()
        return self._adapter_inputs

    def get_inputs(self, k, tau=None, input_selection=None):
        inputs = self.processor.process(
            self.get_adapter_inputs(),
            k, tau,
            selection=input_selection
        )
        return inputs

    def get_outputs(self, k, tau=None, source_function_selection=None):
        outputs = TargetPostProcessor().process(
            self.adapter.get_outputs(source_function_selection),
            k, tau
        )
        return outputs


class SourceFileDataProvider(AdapterProvider):
    def __init__(self, source_file, k_standard=None, compat_mode=False):
        assert not compat_mode

        super().__init__(adapter=SourceFileAdapter(source_file))


class CLASSDataProvider(AdapterProvider):
    def __init__(self, cosmo):
        super().__init__(adapter=CLASSAdapter(cosmo))

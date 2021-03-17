import os

import classynet.pipeline


def get_input_transformer():
    return pipeline.CanonicalInputTransformer()

def get_input_normalizer(normalization_file):
    return pipeline.CanonicalInputNormalizer.from_path(normalization_file)


def get_target_transformer(k):
    return pipeline.CanonicalTargetTransformer(k=k)

def get_target_normalizer_abs_max(normalization_file):
    return pipeline.AbsMaxNormalizer.from_path(
            normalization_file,
            )

def get_target_normalizer(normalization_file):
    return get_target_normalizer_abs_max(normalization_file)


def get_input_transformer_normalizer(normalization_file):
    transformer = get_input_transformer()
    normalizer = get_input_normalizer(normalization_file)
    return pipeline.compose(transformer, normalizer)

def get_target_transformer_normalizer(normalization_file, k):
    transformer = get_target_transformer(k)
    normalizer = get_target_normalizer(normalization_file)
    return pipeline.compose(transformer, normalizer)

def get_pair(normalization_file, k):
    """Return a tuple of (input transformer+normalizer), (output transformer+normalizer)."""
    return get_input_transformer_normalizer(normalization_file), get_target_transformer_normalizer(normalization_file, k)

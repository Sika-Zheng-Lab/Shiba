"""Lightweight regression option mapping for pipeline/container launchers."""

def config_arguments(config):
    """Single option mapping shared by the Python and Snakemake launchers."""
    if config.get("stat_method", "legacy") != "beta-binomial":
        return []
    args = ["--stat-method", "beta-binomial"]
    for key in ("sample_metadata", "formula", "contrast_file", "prediction_grid", "effect_name",
                "min_effect", "p_adjust", "rscript", "event_batch_size", "min_sample_fraction"):
        if config.get(key) is not None:
            args.extend(["--" + key.replace("_", "-"), str(config[key])])
    for key, option in (("coef", "--coef"), ("categorical", "--categorical"), ("continuous", "--continuous"),
                        ("filter_group", "--filter-group")):
        values = config.get(key, [])
        if isinstance(values, str):
            values = [values]
        for value in values:
            args.extend([option, str(value)])
    for factor, level in config.get("reference_levels", {}).items():
        args.extend(["--reference-level", f"{factor}={level}"])
    return args

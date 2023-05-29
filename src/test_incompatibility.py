from loguru import logger


def incompatbility(combined_depth_mash_df):
    """Quick heuristic check whether there is likely differences between long and short read sets

    :param combined_depth_mash_df: combined output dataframe
    :logger: logger
    :return:
    """

    # count
    count = combined_depth_mash_df.loc[
        (combined_depth_mash_df["contig"] != "chromosome")
        & (combined_depth_mash_df["circularity"] == "not_circular")
        & (combined_depth_mash_df["PLSDB_hit"] == "")
    ].shape[0]
    if count >= 5:
        message = (
            "WARNING: "
            + str(count)
            + " non-circular contigs with no PLSDB mash hits were detected. \nThis indicates your long and short read sets may come from different bacterial isolates. \nPlease check this!"
        )
        logger.info(message)

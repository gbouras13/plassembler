import log

def incompatbility(combined_depth_mash_df, logger):
    """Quick heuristic check whether there is likely differences between long and short read sets

    :param combined_depth_mash_df: combined output dataframe
    :logger: logger
    :return:
    """

    # not circular
    combined_depth_mash_df['contig'] != 'chromosome'
    combined_depth_mash_df['circularity'] != 'circular'
    combined_depth_mash_df['PLSDB_hit'] != 'No_hit'

    # count
    count = combined_depth_mash_df.loc[(combined_depth_mash_df['contig'] != 'chromosome') & (combined_depth_mash_df['circularity'] != 'circular') & (combined_depth_mash_df['PLSDB_hit'] != 'No_hit') ].shape[0]

    if count >= 5:
        message = str(count) + ' non-contigs with no PLSDB mash hits were detected. \nThis indicates your long and short read sets may come from different bacterial isolates. \nPlease check this!'
        log.write_message(message, logger)    

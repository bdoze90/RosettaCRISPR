"""This file stores all the information on the indeces that need to be iterated over and structures directing which
structures need which calculations.  Useful as import into every file used in the project."""

OFF_COMBOS_HSU_SPCAS9 = {1957: (1900, 1956, 0, 0),
                           2015: (1958, 2014, 0, 0),
                           1899: (1842, 1898, 0, 0),
                           2073: (2016, 2072, 0, 0),
                           2074: (2075, 2084, 0, 0),
                           2090: (2091, 2102, 0, 0),
                           2103: (2104, 2136, 0, 0),
                           2137: (2138, 2143, 2252, 2257),
                           2144: (2145, 2162, 0, 0),
                           2163: (2164, 2202, 0, 0),
                           2203: (2204, 2209, 2258, 2263),
                           2210: (2211, 2222, 0, 0),
                           2223: (2224, 2241, 2264, 2272),
                           2242: (2243, 2251, 2273, 2281),
                           2282: (2283, 2297, 0, 0)
                           }

OFFBASES_HSU_SPCAS9 = [1899, 1957, 2015, 2073, 2074, 2090, 2103, 2137, 2144, 2163, 2203, 2210, 2223, 2242, 2282]

OFF_COMBOS_KLEINSTIVER = list()

CS_DICT = { "4UN3": {"ChainA": ('r', 0, 81, '', ''),
                     "ChainB": ('protein', '', '', '', ''),
                     "ChainC": ('d', 0, 'n', 'rc', 'TGGTATTG'),
                     "ChainD": ('d', 17, 'n', '', 'TGGTATTG')},

            "4UN4": {"ChainA": ('r', 0, 81, '', ''),
                     "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                     "ChainC": ('d', 17, 'n', 'rc', 'TGGTATTG'),
                     "ChainD": ('d', 18, 'n', '', 'TGGTATTG'),
                     "ChainE": ('d', 0, 17, 'rc', '')},

            # Crystal structure has lost a portion of the sequence and therefore only good for short guides
            "4UN5": {"ChainA": ('r', 0, 81, '', ''),
                     "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                     "ChainC": ('d', 20, 'n', 'rc', 'TGGTATTG'),
                     "ChainD": ('d', 18, 'n', '', 'TGGTATTG'),
                     "ChainE": ('d', 0, 17, 'rc', '')},

            "5FQ5": {"ChainA": ('r', 0, 81, '', ''),
                     "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                     "ChainC": ('d', 19, 'n', 'rc', 'TGGTATTG'),
                     "ChainD": ('d', 18, 'n', '', 'TGGTATTG'),
                     "ChainE": ('d', 0, 17, 'rc', '')},

            # Needs rna_seqs2.txt b/c of different scaffold RNA
            "5F9R": {"ChainA": ('r', 1, 116, '', ''),
                     "ChainB": ('protein', 'n','n','n','n'),
                     "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                     "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},


            # Beginning of the saCas9 structures
            "5CZZ": {"ChainA": ('r', 0, 'n', '', ''),
                     "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                     "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                     "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

            "5AXW": {"ChainA": ('r', 0, 'n', '', ''),
                     "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                     "ChainC": ('d', 0, 'n', 'rc', 'TTGGGTAG'),
                     "ChainD": ('d', 'n', 'n', '', 'TTGGGTAG')},


            # Beginning of the Cas12 structures
            "5XUU": {"ChainA": ('r', 1, 116, '', ''),
                     "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                     "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                     "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

            "5XUS": {"ChainA": ('r', 1, 116, '', ''),
                     "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                     "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                     "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

            "5XUT": {"ChainA": ('r', 1, 116, '', ''),
                     "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                     "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                     "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

            "5XH6": {"ChainA": ('r', 1, 116, '', ''),
                     "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                     "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                     "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

            "5XH7": {"ChainA": ('r', 1, 116, '', ''),
                     "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                     "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                     "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')}
            }

# THIS IS THE ONLY THING THAT NEEDS TO BE CHANGED WHEN SWITCHING COMPUTERS!!!
# IDs:
# MacBookPro 2013: 1
# Linux Lab Computer: 2
# MacBookPro 2017: 3
# BRC Mac Pro: 4
COMPUTER_ID = 4
import logging


def setup(debug_level):
    """Setups a global logger with the name 'logger'. 
    This logger can accessed in any function by "log = logging.getLogger('logger')".
    Three different logging levels:
        0 : print no logging messages
        1 : print information messages (i.e. print all messages invoked with "log.info(message)")
        2 : print debug + information messages (i.e. print all messages invoked with "log.info(message)" or "log.debug(message)")
    """ 

    if debug_level == 0   : # no logging
        formatter = logging.Formatter()  
    elif debug_level == 1 : # info 
        formatter = logging.Formatter(fmt='%(message)s')
    elif debug_level == 2 : # debug
        formatter = logging.Formatter(fmt='%(levelname)s - %(module)s: %(message)s')

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger = logging.getLogger('logger')
    if debug_level == 0 :   # no logging
        logger.setLevel(0)
    elif debug_level == 1 : # info
        logger.setLevel(logging.INFO)
    elif debug_level == 2 : # debug
        logger.setLevel(logging.DEBUG)

    logger.addHandler(handler)
    return logger

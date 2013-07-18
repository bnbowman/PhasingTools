import os, logging, subprocess

log = logging.getLogger()

def run_amplicon_assembly(query, reference, args):
    command_args = create__command(query, reference, args)
    log_command( command_args )
    execute_command( command_args )

def create_blasr_command(query, reference, args):
    log.info("Converting supplied arguments into a Blasr commandline")
    command_args = ['blasr', query, reference]
    for arg, value in args.iteritems():
        arg = '-' + str(arg)
        if value is True:
            command_args += [arg]
        else:
            command_args += [arg, str(value)]
    return command_args

def log_command( args ):
    args = list( args )
    command = args.pop(0).capitalize()
    log.info('Executing "%s" with the following options:' % command)

    # Log any positional arguments
    while not args[0].startswith('-'):
        arg = args.pop(0)
        log.info('\tArgument: %s' % arg)

    # Log any remaining options
    while True:
        if len(args) == 0:
            break
        option_parts = [args.pop(0)]
        while len(args) and not args[0].startswith('-'):
            part = args.pop(0)
            option_parts.append( part )
        log.info('\tOption: "%s"' % ' '.join(option_parts))

def execute_command( command_args ):
    command = command_args[0].capitalize()
    log.info('Executing "%s" command as subprocess' % command)
    with open('/dev/null', 'w') as handle:
        subprocess.check_call( command_args, 
                               stdout=handle, 
                               stderr=subprocess.STDOUT )
    log.info("Subprocess finished successfully")

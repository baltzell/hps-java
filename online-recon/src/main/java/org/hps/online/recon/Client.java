package org.hps.online.recon;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.hps.online.recon.ClientCommand.CleanupCommand;
import org.hps.online.recon.ClientCommand.ConfigCommand;
import org.hps.online.recon.ClientCommand.CreateCommand;
import org.hps.online.recon.ClientCommand.ListCommand;
import org.hps.online.recon.ClientCommand.RemoveCommand;
import org.hps.online.recon.ClientCommand.SettingsCommand;
import org.hps.online.recon.ClientCommand.StartCommand;
import org.hps.online.recon.ClientCommand.StopCommand;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 * Client for interacting with the online reconstruction server.
 */
public final class Client {

    /**
     * Package logger.
     */
    private static Logger LOGGER = Logger.getLogger(Client.class.getPackageName());
    
    /**
     * Hostname of the server with default.
     */
    private String hostname = "localhost";
    
    /**
     * Port of the server with default.
     */
    private int port = 22222;

    /**
     * Output file for writing server responses.
     * By default it is null, which results in output being written to the console (System.out).
     */
    private File outputFile;
    
    /**
     * Parser for base options.
     */
    private CommandLineParser parser = new DefaultParser();
       
    /**
     * Map of names to client commands.
     */
    private Map<String, ClientCommand> commandMap = new LinkedHashMap<String, ClientCommand>();

    /**
     * The set of valid client commands.
     */
    private Set<String> commands;
    
    /**
     * The base options (commands have their own Options objects).
     */
    private static Options OPTIONS = new Options();
    static {
        OPTIONS.addOption(new Option("", "help", false, "print help"));
        OPTIONS.addOption(new Option("p", "port", true, "server port"));
        OPTIONS.addOption(new Option("h", "host", true, "server hostname"));
        OPTIONS.addOption(new Option("o", "output", true, "output file (default writes server responses to System.out)"));
    }
   
    /**
     * Class constructor.
     */
    Client() { 
        // Build the command map.
        buildCommandMap();
    }
    
    /**
     * Print the base command usage.
     */
    private void printUsage() {
        final HelpFormatter help = new HelpFormatter();
        final String commands = String.join(" ", commandMap.keySet());
        help.printHelp("Client [options] [command] [command_options]", "Send commands to the online reconstruction server",
                OPTIONS, "Commands: " + commands + '\n'
                    + "Use 'Client [command] --help' for help with a specific command.");
    }
    
    /**
     * Build a map of names to client commands.
     */
    private void buildCommandMap() {
        
        // Add commands.
        addCommand(new CreateCommand());
        addCommand(new StartCommand());
        addCommand(new StopCommand());
        addCommand(new RemoveCommand());
        addCommand(new ListCommand());
        addCommand(new ConfigCommand());
        addCommand(new CleanupCommand());
        addCommand(new SettingsCommand());
        
        // Define set of valid command names.
        commands = commandMap.keySet();
    }
    
    /**
     * Add a client command.
     * @param command The client command to add
     */
    private void addCommand(ClientCommand command) {
        String name = command.getName();
        commandMap.put(name, command);
    }
    
    /**
     * Get a client command by its name.
     * @param name The name of the client command
     * @return The client command or null if does not exist
     */
    private ClientCommand getCommand(String name) {
        return commandMap.get(name);
    }
    
    /**
     * Scan for a valid command on the command line.
     * @param args The arg array
     * @return The position of the command or -1 if not found
     */
    private int scanForCommand(String args[]) {
        int commandPos = -1;
        for (int i = 0; i < args.length; i++) {
            if (commands.contains(args[i])) {
                commandPos = i;
                break;
            }
        }        
        return commandPos;
    }

    /**
     * Run the client using command line arguments
     * @param args The command line arguments
     */
    void run(String args[]) {
        
        // Scan for command.
        // FIXME: This doesn't work properly if an option value is the same as a command name.
        int commandPos = scanForCommand(args);
        
        // No arguments or no valid command was provided.
        if (args.length == 0 || commandPos == -1) {
            // Print usage and exit.
            printUsage();
            System.exit(0);
        }

        // Get name of command from arguments.
        final String commandName = args[commandPos];

        // Create argument array for base client command.
        String baseArgs[] = new String[commandPos];
        System.arraycopy(args, 0, baseArgs, 0, commandPos);
        //System.out.println(Arrays.asList("baseArgs: " + Arrays.asList(baseArgs)));
                
        // Parse base options.
        CommandLine cl;
        try {
            cl = this.parser.parse(OPTIONS, baseArgs, true);
        } catch (ParseException e) {
            throw new RuntimeException("Error parsing arguments", e);
        }
                        
        if (cl.hasOption("p")) {
            this.port = Integer.parseInt(cl.getOptionValue("p"));
            LOGGER.config("Port: " + this.port);
        }
        
        if (cl.hasOption("h")) {
            this.hostname = cl.getOptionValue("h");
            LOGGER.config("Hostname: " + this.hostname);
        }
        
        if (cl.hasOption("o")) {
            this.outputFile = new File(cl.getOptionValue("o"));
            LOGGER.config("Output file: " + this.outputFile.getPath());
        }
        
        // Create arg array for command.
        ClientCommand command = getCommand(commandName);
        int cmdArrayLen = args.length - baseArgs.length - 1;
        String cmdArgs[] = new String[cmdArrayLen];
        System.arraycopy(args, commandPos + 1, cmdArgs, 0, cmdArrayLen);
        //System.out.println("cmdArgs: " + Arrays.asList(cmdArgs));
                
        // Parse command options.
        DefaultParser commandParser = new DefaultParser();
        CommandLine cmdResult = null;
        try {
            cmdResult = commandParser.parse(command.getOptions(), cmdArgs);
        } catch (ParseException e) {
            command.printUsage();
            throw new RuntimeException("Error parsing command options", e);            
        }
        
        // Setup the command parameters from the parsed options.
        command.process(cmdResult);
        
        // Send the command to server.
        LOGGER.info("Sending command " + command.toString());
        send(command);
    }
    
    /**
     * Send a command to the online reconstruction server.
     * @param command The client command to send
     */
    private void send(ClientCommand command) {
        try (Socket socket = new Socket(hostname, port)) {
            // Send command to the server.           
            PrintWriter writer = new PrintWriter(socket.getOutputStream());
            writer.write(command.toString() + '\n');            
            writer.flush();

            // Get server response.
            InputStream is = socket.getInputStream();
            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            String resp = br.readLine();
            
            // Print server response or write to output file.
            PrintWriter pw = null;
            if (this.outputFile != null) {
                pw = new PrintWriter(this.outputFile);
            }
            if (resp.startsWith("{")) {
                // Handle JSON object.
                JSONObject jo = new JSONObject(resp);
                if (pw != null) {
                    pw.write(jo.toString(4));
                } else {
                    System.out.println(jo.toString(4));
                }
            } else if (resp.startsWith("[")) {
                // Handle JSON array.
                JSONArray ja = new JSONArray(resp);
                if (pw != null) {
                    pw.write(ja.toString(4));
                } else {
                    System.out.println(ja.toString(4));
                }
            } else {
                // Response from server isn't valid JSON.
                throw new RuntimeException("Invalid server response: " + resp.toString());
            }
            if (pw != null) {
                LOGGER.info("Wrote server response to: " + this.outputFile.getPath());
                pw.flush();
                pw.close();
            }
        } catch (Exception e) {
            throw new RuntimeException("Client error", e);
        }
    }
    
    /**
     * Run the client from the command line.
     * @param args The argument array
     */
    public static void main(String[] args) {
        Client client = new Client();
        client.run(args);
    }
}

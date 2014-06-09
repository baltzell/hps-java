package org.hps.monitoring;

/**
 * These strings are used to identify ActionEvents in the MonitoringApplication.
 * 
 * @author Jeremy McCormick <jeremym@slac.stanford.edu>
 */
final class MonitoringCommands {

    private MonitoringCommands() {}

    static final String CONNECT = "connect";
    static final String DISCONNECT = "disconnect";
    static final String SAVE_CONNECTION = "saveConnection";
    static final String LOAD_CONNECTION = "loadConnection";
    static final String RESET_CONNECTION_SETTINGS = "resetConnectionSettings";
    static final String RESET_EVENTS = "resetEvents";
    static final String SAVE_PLOTS = "savePlots";
    static final String RESET_DRIVERS = "resetDrivers";
    static final String SET_EVENT_BUILDER = "setEventBuilder";
    static final String SET_EVENT_REFRESH = "setEventRefresh";
    static final String EXIT = "exit";
    static final String LOG_TO_FILE = "logToFile";
    static final String LOG_TO_TERMINAL = "logToTerminal";
    static final String SCREENSHOT = "screenshot";
    static final String EDIT_EVENT_REFRESH = "editEventRefresh";
    static final String UPDATE_TIME = "updateTime";    
    static final String SET_MAX_EVENTS = "setMaxEvents";
    static final String SAVE_LOG_TABLE = "saveLogTable";
    static final String CLEAR_LOG_TABLE = "clearLogTable";    
    static final String PAUSE = "pause";
    static final String RESUME = "resume";
    static final String NEXT = "next";
    static final String SET_LOG_LEVEL = "setLogLevel";
    static final String AIDA_AUTO_SAVE = "aidaAutoSave";
    static final String SAVE_JOB_SETTINGS = "saveJobSettings";
    static final String LOAD_JOB_SETTINGS = "loadJobSettings";
    static final String RESET_JOB_SETTINGS = "resetJobSettings";
    static final String SET_STEERING_RESOURCE = "setSteeringResource";
    static final String SET_STEERING_FILE = "setSteeringFile";
}
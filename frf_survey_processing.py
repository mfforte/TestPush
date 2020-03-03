import uuid
import arcpy
import sys
import os
import datetime
import urllib
import urllib2
from xml.dom import minidom
import numpy as np
import gzip
import timeit
import requests
from xml.etree import ElementTree as et
import math
import logging
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

arcpy.SetProduct('arcinfo')
arcpy.CheckOutExtension('3D')

# log file settings, writes to file and prints to stderr
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

fh = logging.FileHandler(r"\\coe-samgsp01mob\gis\Work\_C045\Scripts\survey_data_scheduled_task\logs\frf_survey_processing.log")
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)


def download_data_tds(tds_url, out_folder):
    """
    Downloads FRF survey data from the CHL THREDDS server, skipping files that already exist in out_folder
    https://chlthredds.erdc.dren.mil/thredds/catalog/frf/survey_temp/catalog.xml
    :param tds_url: chl thredds server url
    :param out_folder: path to save files
    :return: path to downloaded files
    """
    try:
        if not test_url(tds_url):
            logging.error("ERROR: {0} - THREDDS Server Down".format(tds_url))
            logging.info("="*50)
            sys.exit()
        # Path to THREDDS server with catalog.xml
        xmlFile = urllib.urlopen(tds_url).read()
        xmldoc = minidom.parseString(xmlFile)

        # Get all the files available in this directory
        datasets = xmldoc.getElementsByTagName("dataset")
        existing_files = os.listdir("{0}".format(out_folder))

        # loop through each dataset in the xml file.  Reading its name, check if that file has
        # already been downloaded to the designated location.
        for dataset in datasets:
            if dataset.hasAttribute("name"):
                filename = dataset.getAttribute("name")
                if filename == "survey_temp":
                    continue
                dataset_name = dataset.getAttribute("name")
                download_url = 'https://chlthredds.erdc.dren.mil/thredds/fileServer/frf/survey_temp/{0}'.format(
                        dataset.getAttribute("name"))
                filepath = os.path.join(out_folder, dataset_name)
                searchCounter = 0
                for offlineFileName in existing_files:
                    if dataset_name == offlineFileName:
                        searchCounter += 1
                        logging.info("{0} already downloaded.".format(dataset_name))
                if searchCounter == 0:
                    logging.info('INFO: {0} - downloading file...'.format(dataset_name))
                    urllib.urlretrieve(download_url, filepath)
        return out_folder
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)


def process_data(folder):
    """
    Processes data downloaded from the CHL THREDDS server, skipping files that have been previously processed.
    This is the main function that calls the rest of the supporting functions
    :param folder: folder of survey files to process
    :return:
    """
    try:
        global metadata_counter
        global survey_counter
        metadata_counter = 0
        survey_counter = 0
        for filename in os.listdir(folder):
            if filename.split(".")[-1] == 'xml':
                sds_metadata_id = check_metadata_import(metadata_archive_table, filename)
                if sds_metadata_id is False:
                    sds_metadata_id = upload_metadata(os.path.join(folder, filename), "FRF_Metadata", "FDIF2015!")
                    logging.info('SUCCESS: {0} - uploaded metadata file to Metadata Manager'.format(filename))
                    logging.info("-"*50)
                    metadata_counter += 1
        for filename in os.listdir(folder):
            if filename.split(".")[-1] == 'csv':
                logging.info('INFO: processing file: {0}'.format(filename))
                # filetype = "PROFILE"
                data = load_text(os.path.join(folder, filename))
                if not validate_data(data):
                    data = load_text(os.path.join(folder, filename), skip_header=1)
                surveyjobidfk, survey_date, locality_code = file_details(data, filename)
                if check_import_archive(import_archive_table, surveyjobidfk, filename):
                    continue
                replace = replace_survey_check(surveyjobidfk, filename, 'profile')
                if replace == "duplicate":
                    continue
                elif replace is True:
                    delete_old_survey_profile(surveyjobidfk)
                elif replace is False:
                    continue
                else:
                    pass
                insert_survey_points(survey_point_table, data, sds_metadata_id, surveyjobidfk)
                process_survey_job(survey_job_table, survey_point_table, """SURVEYJOBIDFK = '{0}'""".format(surveyjobidfk),
                                   surveyjobidfk, survey_date, sds_metadata_id)
                process_import_archive(import_archive_table, os.path.join(folder, filename),
                                       sds_metadata_id, surveyjobidfk, survey_date, locality_code)
                if check_metadata_to_surveys(metadata_to_surveys_table, sds_metadata_id, surveyjobidfk):
                    logging.info("SUCCESS: {0} successfully loaded into FDIF".format(filename))
                    logging.info("-"*50)
                    survey_counter += 1
                    continue
                process_metadata_to_surveys(metadata_to_surveys_table, sds_metadata_id, surveyjobidfk)
                logging.info("SUCCESS: {0} successfully loaded into FDIF".format(filename))
                logging.info("-"*50)
                survey_counter += 1
            if filename.split(".")[-1] == 'txt':
                logging.info('INFO: processing file: {0}'.format(filename))
                # filetype = "GRID"
                surveyjobidfk = "FRF_{0}_{1}".format(filename.split("_")[2], filename.split("_")[3])
                if check_import_archive(import_archive_table, surveyjobidfk, filename):
                    continue
                replace = replace_survey_check(surveyjobidfk, filename, 'grid')
                if replace == "duplicate":
                    continue
                elif replace is True:
                    delete_old_survey_grid(surveyjobidfk, filename)
                elif replace is False:
                    continue
                else:
                    pass
                locality_code = str(filename.split("_")[3])
                survey_date = datetime.datetime.strptime(str(filename.split("_")[1]), '%Y%m%d')
                logging.info("INFO: SurveyJobIDFK: {0}".format(surveyjobidfk))
                logging.info("INFO: SurveyDate: {0}".format(survey_date))
                raster_name = "C045_S" + str(filename.split("_")[2]) + "_" + str(filename.split("_")[3]) + "_SPNCM_GRD_" + str(filename.split("_")[1])
                data = load_grid_text(os.path.join(folder, filename))
                try:
                    grid_points = create_scratch_grid_points(data)
                except:
                    continue
                tin = create_tin(grid_points, arcpy.SpatialReference(32119))
                frf_raster = tin_to_raster(tin, raster_name)
                contours = raster_to_contours(frf_raster, 1)
                if insert_contours(elevation_contours, contours, surveyjobidfk, sds_metadata_id):
                    insert_raster(survey_grid_archive, frf_raster)
                    process_import_archive_grid(import_archive_table, os.path.join(folder, filename), sds_metadata_id,
                                                surveyjobidfk, survey_date, locality_code)
                    add_raster_to_mosaic(mosaic_dataset, survey_grid_archive, raster_name)
                    update_mosaic_attributes(mosaic_dataset, raster_name, surveyjobidfk, sds_metadata_id, survey_date,
                                             "GRID", filename)
                    if check_metadata_to_surveys(metadata_to_surveys_table, sds_metadata_id, surveyjobidfk):
                        logging.info("SUCCESS: {0} successfully loaded into FDIF".format(filename))
                        logging.info("-"*50)
                        survey_counter += 1
                        continue
                    process_metadata_to_surveys(metadata_to_surveys_table, sds_metadata_id, surveyjobidfk)
                    logging.info("SUCCESS: {0} successfully loaded into FDIF".format(filename))
                    logging.info("-"*50)
                    survey_counter += 1
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)


def replace_survey_check(surveyjobidfk, filename, type):
    """
    Checks file version, if older survey exists, return True
    :param surveyjobidfk:
    :param filename:
    :return:
    """
    try:
        counter = 0
        if type == 'profile':
            version = int(filename.split(".")[0].split("_")[-1][1:9])
            expression = """surveyNumber = '{0}' AND SURVEYTYPE <> 'GRID'""".format(surveyjobidfk)
        if type == 'grid':
            version = int(filename.split(".")[0].split("_")[-3][1:9])
            expression = """surveyNumber = '{0}' AND SURVEYTYPE = 'GRID'""".format(surveyjobidfk)
        with arcpy.da.SearchCursor(import_archive_table, ["surveyNumber", "VERSION"],
                                   where_clause=expression) as cursor:
            for row in cursor:
                old_version = int(row[1][1:9])
                counter += 1
        if counter > 1:
            logging.info('INFO: {0} - Duplicate surveyNumber, returned more than 1 record: {0}'.format(surveyjobidfk))
            logging.info('-'*50)
            return 'duplicate'
        if counter == 1:
            if old_version < version:
                return True
            if old_version > version:
                return False
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)


def delete_old_survey_profile(surveyjobidfk):
    """
    Removes a survey from all FDIF tables.
    :param surveyjobidfk:
    :param filename:
    :return:
    """
    arcpy.env.overwriteOutput = True

    try:
        archive = arcpy.MakeTableView_management(import_archive_table, "archive")
        arcpy.SelectLayerByAttribute_management(archive, "NEW_SELECTION", """SURVEYNUMBER = '{0}' AND SURVEYTYPE = 'PROFILE'""".format(surveyjobidfk))
        arcpy.DeleteRows_management(archive)

        points = arcpy.MakeFeatureLayer_management(survey_point_table, "points")
        arcpy.SelectLayerByAttribute_management(points, "NEW_SELECTION", """SURVEYJOBIDFK = '{0}'""".format(surveyjobidfk))
        arcpy.DeleteFeatures_management(points)

        arcpy.MakeFeatureLayer_management(survey_job_table, "surveyjob")
        arcpy.SelectLayerByAttribute_management(in_layer_or_view="surveyjob", selection_type="NEW_SELECTION",
                                                where_clause="""SURVEYJOBIDPK = '{0}'""".format(surveyjobidfk))
        arcpy.DeleteFeatures_management("surveyjob")

        metadatatosurveys = arcpy.MakeTableView_management(metadata_to_surveys_table, "metadatatosurveys")
        arcpy.SelectLayerByAttribute_management(metadatatosurveys, "NEW_SELECTION", """SURVEYJOBIDFK = '{0}'""".format(surveyjobidfk))
        arcpy.DeleteRows_management(metadatatosurveys)
        logging.info("SUCCESS: deleted survey profile")
        return True
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def delete_old_survey_grid(surveyjobidfk, filename):
    """
    Removes a GRID survey from all FDIF tables.
    :param surveyjobidfk:
    :param filename:
    :return:
    """
    arcpy.env.overwriteOutput = True

    try:
        arcpy.MakeTableView_management(import_archive_table, "archive")
        arcpy.SelectLayerByAttribute_management(in_layer_or_view="archive", selection_type="NEW_SELECTION",
                                                where_clause="""SURVEYNUMBER = '{0}' AND SURVEYTYPE = 'GRID'""".format(surveyjobidfk))
        arcpy.DeleteRows_management("archive")

        arcpy.MakeFeatureLayer_management(elevation_contours, "elevationcontours")
        arcpy.SelectLayerByAttribute_management(in_layer_or_view="elevationcontours", selection_type="NEW_SELECTION",
                                                where_clause="""SURVEYJOBIDFK = '{0}'""".format(surveyjobidfk))
        arcpy.DeleteFeatures_management("elevationcontours")

        arcpy.MakeTableView_management(metadata_to_surveys_table, "metadatatosurveys")
        arcpy.SelectLayerByAttribute_management(in_layer_or_view="metadatatosurveys", selection_type="NEW_SELECTION",
                                                where_clause="""SURVEYJOBIDFK = '{0}'""".format(surveyjobidfk))
        arcpy.DeleteRows_management("metadatatosurveys")

        arcpy.MakeMosaicLayer_management(mosaic_dataset, "mosaic")
        arcpy.RemoveRastersFromMosaicDataset_management("mosaic", """SURVEYJOBIDFK = '{0}'""".format(surveyjobidfk))

        raster_name = "C045_S" + str(filename.split("_")[2]) + "_" + str(filename.split("_")[3]) + "_SPNCM_GRD_" + str(
            filename.split("_")[1])
        arcpy.MakeTableView_management(survey_grid_archive, "surveygrid")
        arcpy.SelectLayerByAttribute_management(in_layer_or_view="surveygrid", selection_type="NEW_SELECTION",
                                                where_clause="""NAME = '{0}'""".format(raster_name))
        arcpy.DeleteRows_management("surveygrid")
        logging.info("SUCCESS: deleted survey grid")
        return True
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def test_url(url):
    """
    Test if url is online, 3 attempts are made to reach url.
    :param url: url to test
    :return: boolean, True if online
    """
    attempts = 0
    success = False
    while attempts < 3 and not success:
        try:
            urllib2.urlopen(url)
            success = True
            return success
        except:
            attempts += 1
            return success


def load_text(filepath, skip_header=0):
    """
    Loads FRF survey text data into a structured numpy array
    :param filepath: path to text survey text file
    :param skip_header: number of lines to skip in text file
    :return: numpy array
    """
    try:
        fields = ['localityField', 'profileField', 'surveyNumField', 'yField', 'xField', 'eastingField', 'northingField',
                  'distOffshoreField', 'distLongshoreField', 'zField', 'ellipsField', 'collectionDateField_temp',
                  'collectionTimeField', 'spmField']
        data = np.genfromtxt(filepath, dtype=None, names=fields, delimiter=',', skip_header=skip_header)
        logging.info("SUCCESS: load_text")
        return data
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)


def load_grid_text(filepath):
    """
    Load txt GRID file into a structured numpy array. Assumes dataType: f8 for fields.
    :param filepath: path to txt file
    :return: numpy array
    """
    try:
        fields = ['x', 'y', 'z']
        data = np.genfromtxt(filepath, names=fields, delimiter=',')
        logging.info("SUCCESS: load_grid_text")
        return data
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        error_msg = exc_type, fname, exc_tb.tb_lineno, e
        logging.exception(error_msg)
        email_alert_error(error_msg)
        sys.exit()


def create_scratch_grid_points(data):
    """
    Creates a temporary point feature class from gridded points.
    Output Spatial Reference: 32119 (North Carolina State Plane Meters)
    :param data: return of load_data function, assumes WGS84 projection
    :return: path to in_memory points feature class
    """
    arcpy.env.overwriteOutput = True
    scratch_table = "grid_points"
    txt_sr = arcpy.SpatialReference(4326)
    fc_sr = arcpy.SpatialReference(32119)
    try:
        arcpy.CreateFeatureclass_management("in_memory", scratch_table, "POINT", has_z="ENABLED", spatial_reference=fc_sr)
        with arcpy.da.InsertCursor("in_memory\\{0}".format(scratch_table), ["SHAPE@"]) as i_cur:
            for i in data:
                xyz = arcpy.PointGeometry(arcpy.Point(i['x'], i['y'], i['z']), txt_sr)
                xyz.projectAs(fc_sr, "WGS_1984_(ITRF00)_To_NAD_1983")
                i_cur.insertRow([xyz])
        logging.info("SUCCESS: created scratch grid points")
        return "in_memory\\{0}".format(scratch_table)
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        error_msg = exc_type, fname, exc_tb.tb_lineno, e
        logging.exception(error_msg)
        email_alert_error(error_msg)
        # sys.exit()


def create_tin(points, spref):
    """
    Create TIN from 3D point feature class.
    requires 3D Analyst extension
    :param points: path to points feature class
    :param spref: spatial reference object, spatial reference of the output tin
    :return: path to tin
    """
    try:
        arcpy.env.overwriteOutput = True
        scratch_folder = arcpy.env.scratchFolder
        tin = arcpy.CreateTin_3d(os.path.join(scratch_folder, "frf_tin"), spref, points + " " + "Shape.Z Mass_Points <None>")
        logging.info("SUCCESS: created TIN")
        return tin
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def tin_to_raster(tin, raster_name):
    """
    Creates a raster from TIN.
    requires 3D Analyst extension
    :param tin: path to tin
    :param raster_name: output raster name
    :return: path to output raster
    """
    try:
        arcpy.env.overwriteOutput = True
        raster = arcpy.TinRaster_3d(tin, "in_memory\\{0}".format(raster_name), "FLOAT", "LINEAR",
                                    "OBSERVATIONS 250", 1)
        logging.info("SUCCESS: tin to raster, name = {0}".format(raster_name))
        return raster
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def raster_to_contours(raster, contour_interval):
    """
    Create contour lines from raster. WGS84 output coord system.
    :param raster: path to raster file
    :param contour_interval: contour interval, units in coordinate system of raster
    :return: path to output contours
    """
    try:
        arcpy.env.overwriteOutput = True
        arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(4326)
        arcpy.env.geographicTransformations = "WGS_1984_(ITRF00)_To_NAD_1983"
        scratch_db = arcpy.env.scratchGDB
        contour_table = "frf_contours"
        minimum = int(math.ceil(arcpy.Raster(raster).minimum))
        contour = arcpy.Contour_3d(raster, os.path.join(scratch_db, contour_table), contour_interval, minimum)
        arcpy.TrimLine_edit(contour, "500 Meters", "DELETE_SHORT")
        arcpy.ClearEnvironment("outputCoordinateSystem")
        arcpy.ClearEnvironment("geographicTransformations")
        logging.info("SUCCESS: raster to contours")
        return contour
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def insert_contours(target_contour_table, source_contour_table, surveyjobidfk, sds_metadata_id):
    """
    Insert contours in FDIF.
    :param target_contour_table: ElevationContours table
    :param source_contour_table: contours to insert
    :param surveyjobidfk: surveyJobIDFK
    :param sds_metadata_id: sdsMetadataID
    :return: boolean, True if successful
    """
    try:
        scratch_table = "contours"
        fc_sr = arcpy.SpatialReference(4326)
        arcpy.CreateFeatureclass_management("in_memory", scratch_table, "POLYLINE", template=elevation_contours,
                                            spatial_reference=fc_sr)
        with arcpy.da.SearchCursor(source_contour_table, ["SHAPE@", "Contour"]) as s_cur:
            with arcpy.da.InsertCursor(target_contour_table, ["SHAPE@", "surveyJobIDFK", "elevationDatum", "elevationUOM",
                                                              "contourElevation", "elevationContourType", "sdsFeatureName",
                                                              "sdsFeatureDescription", "sdsMetadataID", "sdsID",
                                                              "contourElevationIDPK", "elevationContourIDPK"]) as i_cur:
                for row in s_cur:
                    i_cur.insertRow([row[0], surveyjobidfk, "NAVD88", "M", row[1], "INDEX", row[1],
                                    "Derived contour from FRF survey", sds_metadata_id,
                                    '{' + str(uuid.uuid4()).upper() + '}', '{' + str(uuid.uuid4()).upper() + '}',
                                    '{' + str(uuid.uuid4()).upper() + '}'])
        arcpy.Append_management("in_memory\\{0}".format(scratch_table), elevation_contours, "NO_TEST")
        logging.info("SUCCESS: ElevationContours records inserted")
        return True
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def insert_raster(target_raster_catalog, raster):
    """
    Inserts raster into Raster Catalog
    :param target_raster_catalog: raster catalog table
    :param raster: raster file to insert
    :return:
    """
    try:
        arcpy.RasterToGeodatabase_conversion(raster, target_raster_catalog)
        logging.info("SUCCESS: insert raster in raster catalog")
        return True
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)


def add_raster_to_mosaic(mosaic_dataset, survey_grid_archive, raster_name):
    """
    Add raster to mosaic dataset from Raster Catalog
    :param mosaic_dataset: mosaic dataset table
    :param survey_grid_archive: raster catalog
    :param raster_name: name of raster in raster catalog
    :return:
    """
    try:
        arcpy.AddRastersToMosaicDataset_management(mosaic_dataset, "Table", survey_grid_archive, "NO_CELL_SIZES", "NO_BOUNDARY",
                                                   "NO_OVERVIEWS", "", "0", "1500", "", "Name='" + raster_name + "'",
                                                   "NO_SUBFOLDERS", "EXCLUDE_DUPLICATES", "BUILD_PYRAMIDS",
                                                   "CALCULATE_STATISTICS", "NO_THUMBNAILS", "",
                                                   "NO_FORCE_SPATIAL_REFERENCE")
        logging.info("SUCCESS: {0} - raster addded to mosaic dataset".format(mosaic_dataset))
    except Exception as e:
        logging.exception(e)
        #email_alert_error(e.message)
        sys.exit()


def update_mosaic_attributes(mosaic_dataset, raster_name, surveyjobidfk, sds_metadata_id, survey_date, survey_type, filename):
    """
    Updates the mosaic dataset attributes
    :param mosaic_dataset: mosaic dataset table
    :param raster_name: name of raster in mosaic dataset, found in NAME field
    :param surveyjobidfk: surveyJobIDFK of survey
    :param sds_metadata_id: sdsMetadataID of survey
    :param survey_date: date survey was collected; datatime object
    :param survey_type: type of survey, in filename of survey to process
    :return:
    """
    try:
        expression = """NAME = '{0}'""".format(raster_name)
        with arcpy.da.UpdateCursor(mosaic_dataset, ["SURVEYDATE", "SURVEYJOBIDFK", "SURVEYTYPE", "PRODUCTNAME",
                                                    "SDSMETADATAID", "MAXPS", "MINPS", "NAME"], where_clause=expression) as u_cur:
            for row in u_cur:
                row[0] = survey_date
                row[1] = surveyjobidfk
                row[2] = survey_type
                row[3] = filename
                row[4] = sds_metadata_id
                row[5] = 100
                row[6] = 0
                u_cur.updateRow(row)
        logging.info("SUCCESS: updated mosaic dataset attributes")
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def insert_survey_points(survey_point_table, data, sds_metadata_id, surveyjobidfk):
    """
    Inserts FRF Survey data into the SurveyPoint table
    :param survey_point_table: SurveyPoint table in FDIF
    :param data: return from load_text/load_grid_text function
    :param sds_metadata_id: sdsMetadataID of survey
    :return: surveyjobidfk; False if unsuccessful
    """
    scratch_table = "profile_points"
    fc_sr = arcpy.SpatialReference(4326)
    try:
        arcpy.CreateFeatureclass_management("in_memory", scratch_table, "POINT", template=survey_point_table,
                                            spatial_reference=fc_sr)
        for x in data:
            with arcpy.da.InsertCursor("in_memory\\{0}".format(scratch_table), ["SHAPE@XY", "profileID", "planarCoordSys",
                                                                       "elevationDatum", "sourceType", "elevationUoM",
                                                                       "localityCode", "xLocation", "yLocation",
                                                                       "Easting", "Northing", "surveyPointElev",
                                                                       "ellipsoidHeight", "secondsPastMidnight",
                                                                       "collectionTime", "FRFLongshoreCoord",
                                                                       "FRFOffshoreCoord", "xLocationUoM",
                                                                       "yLocationUoM", "sdsFeatureDescription",
                                                                       "surveyJobIDFK", "collectionDate",
                                                                       "surveyPointIDPK", "sdsFeatureName",
                                                                       "sdsID", "sdsMetadataID"]) as icur:
                icur.insertRow([(x['xField'], x['yField']), x['profileField'], "North Carolina, State Plane, Meter",
                                "NAVD_88", "other", "meter", x['localityField'], x['xField'], x['yField'],
                                x['eastingField'], x['northingField'], x['zField'], x['ellipsField'],
                                x['spmField'], x['collectionTimeField'], x['distLongshoreField'],
                                x['distOffshoreField'], "degrees", "degrees",
                                "FRF Survey - {0}".format(str(x['collectionDateField_temp'])),
                                surveyjobidfk, datetime.datetime.strptime(str(x['collectionDateField_temp']), '%Y%m%d'),
                                datetime.datetime.now().strftime("%Y%m%d%H%M%S%f"), '{0}_{1}_{2}'.format(x['surveyNumField'],
                                                                                                x['profileField'],
                                                                                                x['zField']),
                                '{' + str(uuid.uuid4()).upper() + '}', sds_metadata_id])
        arcpy.Append_management("in_memory\\{0}".format(scratch_table), survey_point_table, "NO_TEST")
        logging.info("SUCCESS: SurveyPoint loaded")
        return True
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        return False


def process_survey_job(survey_job_table, clause_table, where_clause, surveyjobidfk, survey_date, sds_metadata_id):
    """
    Inserts record in the SurveyJob feature class for each survey
    :param survey_job_table: SurveyJob table
    :param clause_table: table that contains the data to compute the minimum bounding geometry tool
    :param where_clause: SQL clause to limit points, etc., applied to clause_table param
    :param surveyjobidfk: surveyJobIDFK
    :param survey_date: datetime of survey, gathered from filename
    :param sds_metadata_id: sdsMetadataID of survey
    :return: boolean
    """
    try:
        scratch_table = "survey_job"
        fc_sr = arcpy.SpatialReference(4326)
        arcpy.CreateFeatureclass_management("in_memory", scratch_table, "POLYGON", template=survey_job_table,
                                            spatial_reference=fc_sr)

        feat_lyr = arcpy.MakeFeatureLayer_management(clause_table, "feat_lyr", where_clause=where_clause)
        boundary = arcpy.MinimumBoundingGeometry_management(feat_lyr, "in_memory\surveyjob", "RECTANGLE_BY_AREA", "ALL")
        with arcpy.da.InsertCursor("in_memory\\{0}".format(scratch_table), ["SHAPE@", "surveyJobIDPK", "surveyAgency",
                                                      "surveyType","surveyAvailability", "sourceDataLocation",
                                                      "sourceDataFormat","surveyDateStart", "surveyDateEnd",
                                                      "sdsMetadataID","sdsFeatureName", "sdsFeatureDescription",
                                                                            "sdsID"]) as icur:
            for row in arcpy.da.SearchCursor(boundary, ["SHAPE@"]):
                icur.insertRow([row[0], surveyjobidfk, "FRF", "controlledSurvey", "complete", "DMZ FRF Schema",
                                "GIS Point Feature Class", survey_date, survey_date, sds_metadata_id,
                                '{0} (FRF)'.format(surveyjobidfk),
                                "Field Research Facility - LARC Survey {0}".format(surveyjobidfk),
                                '{' + str(uuid.uuid4()).upper() + '}'])
        arcpy.Append_management("in_memory\\{0}".format(scratch_table), survey_job_table, "NO_TEST")
        logging.info("SUCCESS: SurveyJob loaded")
        return True
    except Exception as e:
        logging.error("ERROR: SurveyJob load failed: {0}".format(surveyjobidfk))
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def process_import_archive(import_archive_table, filepath, sds_metadata_id, surveyjobidfk, survey_date, thetitle):
    """
    Insert record into the FRFImportArchive table
    :param import_archive_table: FRFImportArchive table
    :param filepath: path to gzip file containing survey data
    :param sds_metadata_id: sdsMetadataID of survey
    :param surveyjobidfk: surveyJobIDFK
    :param survey_date: datetime of survey, gathered from filename
    :return: boolean
    """
    try:
        zipped_file = gzip_file(filepath)
        read_gzip = open(zipped_file, "rb").read()
        filename = os.path.basename(filepath)
        version = filename.split("_")[-1].split(".")[0]
        with arcpy.da.InsertCursor(import_archive_table, ["textFilename", "loadDate", "surveyDate", "surveyNumber",
                                                     "surveyType", "sdsMetadataID", "contourQC", "rasterQC", "pointsQC",
                                                     "boundaryQC", "profileQC", "digitalTextFile", "METHOD", "THETITLE",
                                                          "VERSION"]) as icur:
            icur.insertRow([filename, datetime.datetime.now(), survey_date, surveyjobidfk, "PROFILE", sds_metadata_id,
                            "N/A", "N/A", "NO", "NO", "NO", read_gzip, "GPS", thetitle, version])
        os.remove(zipped_file)
        logging.info("SUCCESS: FRFImportArchive record inserted")
        return True
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def process_import_archive_grid(import_archive_table, filepath, sds_metadata_id, surveyjobidfk, survey_date, thetitle):
    """
    Insert record into the FRFImportArchive table for GRID files (.txt) files
    :param import_archive_table: FRFImportArchive table
    :param filepath: path to gzip file containing survey data
    :param sds_metadata_id: sdsMetadataID of survey
    :param surveyjobidfk: surveyJobIDFK
    :param survey_date: datetime of survey, gathered from filename
    :return: boolean
    """
    try:
        zipped_file = gzip_file(filepath)
        read_gzip = open(zipped_file, "rb").read()
        filename = os.path.basename(filepath)
        version = filename.split("_")[-3]
        with arcpy.da.InsertCursor(import_archive_table, ["textFilename", "loadDate", "surveyDate", "surveyNumber",
                                                     "surveyType", "sdsMetadataID", "contourQC", "rasterQC", "pointsQC",
                                                     "boundaryQC", "profileQC", "digitalTextFile", "METHOD", "THETITLE",
                                                          "VERSION"]) as icur:
            icur.insertRow([filename, datetime.datetime.now(), survey_date, surveyjobidfk, "GRID", sds_metadata_id,
                            "NO", "NO", "N/A", "NO", "N/A", read_gzip, "GPS", thetitle, version])
        os.remove(zipped_file)
        logging.info("SUCCESS: FRFImportArchive record inserted")
        return True
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def process_metadata_to_surveys(metadata_to_surveys_table, sds_metadata_id, surveyjobidfk):
    """
    Inserts record into the FRFMetadataToSurveys table
    :param metadata_to_surveys_table: MetadataToSurveys table in FDIF
    :param sds_metadata_id: sdsMetadataID
    :param surveyjobidfk: surveyJobIDFK
    :return: boolean
    """
    try:
        with arcpy.da.InsertCursor(metadata_to_surveys_table, ["sdsMetadataID", "surveyJobIDFK"]) as icur:
            icur.insertRow([sds_metadata_id, surveyjobidfk])
        logging.info("INFO: FRFMetadataToSurveys record inserted")
        return True
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def check_import_archive(import_archive_table, surveyidfk, text_filename):
    """
    Test if the FRFImportArchive table is populated in FDIF.
    :param import_archive_table: FRFImportArchive table in FDIF
    :param surveyidfk: surveyJobIDFK
    :param text_filename: filename of survey including file extension
    :return: boolean
    """
    try:
        counter = 0
        if surveyidfk.split("_")[-1] == "FRF":
            surveyidfk = "_".join(surveyidfk.split("_")[0:2])
            expression = """textFilename = '{0}' AND surveyNumber LIKE '{1}%'""".format(text_filename, surveyidfk)
        else:
            expression = """textFilename = '{0}' AND surveyNumber = '{1}'""".format(text_filename, surveyidfk)
        with arcpy.da.SearchCursor(import_archive_table, ["textFilename", "surveyNumber"],
                                   where_clause=expression) as cursor:
            for row in cursor:
                counter += 1
        if counter > 0:
            logging.info('INFO: {0} - FRFImportArchive previously loaded, skipping...'.format(text_filename))
            logging.info('-'*50)
            return True
        else:
            return False
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def check_metadata_to_surveys(metadata_to_surveys_table, sds_metadata_id, surveyjobidfk):
    """
    Test if record is inserted in the MetadataToSurveys table in FDIF.
    :param metadata_to_surveys_table: FRFMetadataToSurveys table
    :param surveyidfk: SurveyJobIDFK
    :param sds_metadata_id: sdsMetadataID
    :return: boolean, True if record already exists
    """
    try:
        counter = 0
        expression = """SDSMETADATAID = '{0}' AND SURVEYJOBIDFK = '{1}'""".format(sds_metadata_id, surveyjobidfk)
        with arcpy.da.SearchCursor(metadata_to_surveys_table, ["SDSMETADATAID", "SURVEYJOBIDFK"],
                                   where_clause=expression) as cursor:
            for row in cursor:
                counter += 1
        if counter > 0:
            logging.info("INFO: FRFMetadataToSurveys previously loaded, skipping..")
            return True
        else:
            return False
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def check_metadata_import(metadata_archive_table, filepath):
    """
    Test if metadata file is loaded in FDIF. If it is loaded, make sure it is uploaded to Metadata Manager.
    :param metadata_archive_table: FRFMetadataArchive table
    :param filepath: path to metadata (xml) file
    :return: sdsMetadataID if metadata is already imported, False if it is not loaded
    """
    try:
        filename = os.path.basename(filepath)
        counter = 0
        expression = """origFilename = '{0}'""".format(filename)
        with arcpy.da.SearchCursor(metadata_archive_table, ["origFilename", "sdsMetadataID"],
                                   where_clause=expression) as cursor:
            for row in cursor:
                counter += 1
                sds_metadata_id = row[1]
        if counter > 0:
            logging.info("INFO: {0} - metadata file already loaded in FDIF".format(filename))
            if check_metadata_manager(sds_metadata_id) == sds_metadata_id:
                logging.info("-"*50)
                return sds_metadata_id
            else:
                logging.warning("WARNING: {0} - update sdsMetadataID in the database".format(sds_metadata_id))
                return False
        if counter == 0:
            return False
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def check_metadata_manager(sds_metadata_id):
    """
    Check if metadata file is loaded in Metadata Manager based on the sdsMetadataID
    :param sds_metadata_id: sdsMetadataID of survey
    :return: sdsMetadataID
    """
    metadata_query = 'http://metadata.usace.army.mil/geoportal/rest/find/document?f=xjson&searchText=uuid:"{0}"'.format(sds_metadata_id)
    r = requests.get(metadata_query)
    data = r.json()
    total_results = data['totalResults']
    if total_results == 1:
        uuid = data['features'][0]['properties']['catalog.docuuid']
        sds_metadata_id = str(uuid)
        return sds_metadata_id
    else:
        return False


def upload_metadata(xmlfile, username, password):
    """
    Uploads the metadata file (.xml) to Metadata Manager
    :param xmlfile: path to xml file
    :param username: username to Metadata Manager
    :param password: password to Metadata Manager
    :return: sdsMetadata if successful, False if unsuccessful
    """
    try:
        metadata_manager = 'http://metadata.usace.army.mil/geoportal/rest/manage/document'
        xmlfile = open(xmlfile, 'rb')
        tree = et.parse(xmlfile)
        metadata_title = tree.find(".//{http://www.isotc211.org/2005/gmd}identificationInfo/{http://www.isotc211.org/2005/"
                             "gmd}MD_DataIdentification/{http://www.isotc211.org/2005/gmd}citation/{http://"
                             "www.isotc211.org/2005/gmd}CI_Citation/{http://www.isotc211.org/2005/gmd}title/"
                             "{http://www.isotc211.org/2005/gco}CharacterString").text
        # upload metadata record to Metadata Manager
        r = requests.put(metadata_manager, data=xmlfile, auth=(username, password))
        if r.status_code == requests.codes.created:
            date_uploaded = datetime.datetime.date(datetime.datetime.now())
            sds_metadata_id = get_metadata_uuid(metadata_title, date_uploaded)
            return sds_metadata_id
        elif r.status_code == requests.codes.ok:
            date_uploaded = datetime.datetime.date(datetime.datetime.now())
            sds_metadata_id = get_metadata_uuid(metadata_title, date_uploaded)
            return sds_metadata_id
        elif r.raise_for_status is not None:
            r.raise_for_status
            return False
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        return False


def get_metadata_uuid(metadata_title, date_uploaded):
    """
    Return docuuid of recently uploaded metadata file from Metadata Manager based on title.
    :param metadata_title: title of the metadata file
    :param date_uploaded: datetime the file was uploaded
    :return: sdsMetadataID of the metadata record found
    """
    metadata_query = 'http://metadata.usace.army.mil/geoportal/rest/find/document?f=xjson&searchText=title:"{0}"'.format(metadata_title)
    r = requests.get(metadata_query)
    data = r.json()
    total_results = data['totalResults']
    if total_results == 1:
        uuid = data['features'][0]['properties']['catalog.docuuid']
        sds_metadata_id = str(uuid)
        return sds_metadata_id
    elif total_results > 1:
        for f in data['features']:
            date = datetime.datetime.strptime(f['properties']['catalog.inputdate'][0:10], "%Y-%m-%d")
            if date == date_uploaded:
                sds_metadata_id = f['properties']['catalog.docuuid']
                return sds_metadata_id


def validate_data(data):
    """
    Check the structured array (latitude column) to see if a header was used in the file.
    If a header is detected, re-run load_data module with skip_header = 1
    :param data: out from load_data function
    :return: boolean
    """
    try:
        int(data[0][3])
        return True
    except ValueError:
        return False


def file_details(data, filename):
    """
    Returns details about the file
    :param data:
    :param filename: filename of file being processed
    :return:
    """
    try:
        for x in data:
            locality_code = str(x['localityField'])
            if len(str(x['surveyNumField'])) == 4:
                surveyjobidfk = 'FRF_{0}_{1}'.format(str(x['surveyNumField']), locality_code)
            else:
                # surveyjobidfk = 'FRF_{0}_{1}'.format(str(filename.split("_")[2]), str(filename.split("_")[3]))
                surveyjobidfk = 'FRF_{0}_{1}'.format(str(filename.split("_")[2]), locality_code)
            if len(surveyjobidfk.split("_")[1]) != 4:
                raise Exception('surveyNumber is not of length 4')
            survey_date = datetime.datetime.strptime(str(x['collectionDateField_temp']), '%Y%m%d')
            break
        logging.info("INFO: surveyJobIDFK: {0}".format(surveyjobidfk))
        logging.info("INFO: surveyDate: {0}".format(survey_date))
        return surveyjobidfk, survey_date, locality_code
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def gzip_file(filepath):
    """
    GZIP a file
    :param filepath: path to file
    :return: path to gzip file
    """
    try:
        f_in = open(filepath, 'rb')
        f_out = gzip.open(filepath + '.gz', 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        file_location = filepath + '.gz'
        return file_location
    except Exception as e:
        logging.exception(e)
        email_alert_error(e.message)
        sys.exit()


def email_alert():
    try:
        msg = MIMEMultipart()
        smtp = smtplib.SMTP()
        smtp.set_debuglevel(0)
        smtp.connect('')
        to_email = ["Michael.F.Forte@erdc.dren.mil"]
        from_email = 'NOREPLY@usace.army.mil'
        msg['Subject'] = 'FRF Report'
        msg['From'] = from_email
        date = datetime.datetime.now().strftime("%m/%d/%Y %H:%M")
        text = 'DATE: {0}\nFROM: {1}\nTO: {2}\nSUBJECT: {3}\n\n\nSTATS:\n {4} metadata files processed\n {5} surveys ' \
               'processed'.format(date, from_email, to_email, msg['Subject'], str(metadata_counter), str(survey_counter))
        body = MIMEText(text, 'plain')
        msg.attach(body)
        smtp.sendmail(from_email, to_email, msg.as_string())
        smtp.quit()
        logging.info("email sent successfully")
    except Exception as e:
        logging.exception(e)
        sys.exit()


def email_alert_error(error_message):
    try:
        msg = MIMEMultipart()
        smtp = smtplib.SMTP()
        smtp.set_debuglevel(0)
        smtp.connect('')
        to_email = ["Michael.F.Forte@erdc.dren.mil"]
        from_email = 'NOREPLY@usace.army.mil'
        msg['Subject'] = 'FRF Report'
        msg['From'] = from_email
        date = datetime.datetime.now().strftime("%m/%d/%Y %H:%M")
        text = 'DATE: {0}\nFROM: {1}\nTO: {2}\nSUBJECT: {3}\n\n\nERROR:\n {4}'.format(date, from_email, to_email,
                                                                                      msg['Subject'], error_message)
        body = MIMEText(text, 'plain')
        msg.attach(body)
        smtp.sendmail(from_email, to_email, msg.as_string())
        smtp.quit()
        logging.info("email sent successfully")
        sys.exit()
    except Exception as e:
        logging.exception(e)
        sys.exit()

if __name__ == "__main__":
    tic = timeit.default_timer()
    # Environment settings
    arcpy.SetLogHistory(False)
    arcpy.env.overwriteOutput = True
    arcpy.env.maintainSpatialIndex = True
    arcpy.env.autoCommit = 5000
    # development db
    # arcpy.env.workspace = r"\\coe-samgsp01mob\gis\Work\_C045\Scripts\survey_data_scheduled_task\development_db\FRF.gdb"
    arcpy.env.workspace = r"\\coe-samgsp01mob\GIS\Tools\SDE\administrator\RedDragon_COE_GEO_FRFD.sde"
    survey_point_table = "SurveyPoint"
    survey_job_table = "SurveyJob"
    import_archive_table = "FRFImportArchive"
    metadata_to_surveys_table = "FRFMetadataToSurveys"
    metadata_archive_table = "FRFMetadataArchive"
    elevation_contours = "ElevationContours"
    survey_grid_archive = "SPNCM_SURVEY_GRID_ARCHIVE"
    mosaic_dataset = "WGS84_SURVEY_GRID"
    tds_url = 'https://chlthredds.erdc.dren.mil/thredds/catalog/frf/survey_temp/catalog.xml'
    out_folder = r'\\coe-samgsp01mob\gis\Work\_C045\Scripts\survey_data_scheduled_task\downloaded_files'

    # processing_folder = out_folder
    processing_folder = download_data_tds(tds_url, out_folder)
    process_data(processing_folder)
    toc = timeit.default_timer()
    email_alert()
    logging.info('STATS: elapsed time: {0} seconds'.format(round(toc - tic, 2)))
    logging.info('STATS: {0} metadata files processed'.format(str(metadata_counter)))_al
    logging.info('STATS: {0} surveys processed'.format(str(survey_counter)))
    logging.info("="*50)
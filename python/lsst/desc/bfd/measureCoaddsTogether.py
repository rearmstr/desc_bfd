
from lsst.afw.table import SourceCatalog
from lsst.pipe.base import (
    CmdLineTask, ArgumentParser,
)
from lsst.meas.base import NoiseReplacerConfig, NoiseReplacer
from lsst.pex.config import Config, ConfigField, Field

from lsst.coadd.utils.coaddDataIdContainer import ExistingCoaddDataIdContainer
from lsst.pipe.tasks.multiBand import MergeSourcesRunner, getShortFilterName


__all__ = ("ProcessCoaddsTogetherConfig", "ProcessCoaddsTogetherTask")


class ProcessCoaddsTogetherConfig(Config):
    images = Field(
        doc="Coadd image DatasetType used as input (one for each band)",
        default="deepCoadd_calexp",
        dtype=str,
    )
    ref = Field(
        doc="Coadd catalog DatasetType reference input (one instance across all bands).",
        default="deepCoadd_ref",
        dtype=str,
    )
    output = Field(
        doc="Output catalog DatasetType (one instance across all bands)",
        default=None,   # Must be overridden by derived classes to a DatasetType known to obs_base
        dtype=str,
    )
    deblendReplacer = ConfigField(
        dtype=NoiseReplacerConfig,
        doc=("Details for how to replace neighbors with noise when applying deblender outputs. "
             "Ignored if `useDeblending == False`.")
    )
    deblendCatalog = Field(
        doc=("Catalog DatasetType from which to extract deblended [Heavy]Footprints (one for each band). "
             "Ignored if 'useDeblending == False'."),
        default="deepCoadd_meas",
        dtype=str,
    )


class ProcessCoaddsTogetherTask(CmdLineTask):
    _DefaultName = "processCoaddsTogether"
    ConfigClass = ProcessCoaddsTogetherConfig

    # This feeds the runDataRef() method all bands at once, rather than each one separately.
    # The name reflects how it's used elsewhere, not what it does
    RunnerClass = MergeSourcesRunner

    # TODO: override DatasetType introspection for PipelineTask.  Probably
    # blocked on DM-16275.

    @classmethod
    def _makeArgumentParser(cls):
        # Customize argument parsing for CmdLineTask.
        parser = ArgumentParser(name=cls._DefaultName)
        # This should be config.images.name, but there's no way to pass that
        # information in here in Gen2.
        datasetType = "deepCoadd_calexp"
        parser.add_id_argument("--id", datasetType,
                               ContainerClass=ExistingCoaddDataIdContainer,
                               help="data ID, e.g. --id tract=12345 patch=1,2 filter=g^r^i")
        return parser

    def __init__(self, *, config=None, refSchema=None, butler=None, initInputs=None, **kwds):
        super().__init__(config=config, **kwds)
        if refSchema is None:
            if butler is None:
                if initInputs is not None:
                    refSchema = initInputs.get("refSchema", None)
                if refSchema is None:
                    refSchema = SourceCatalog.Table.makeMinimalSchema()
            else:
                refSchema = butler.get(self.config.ref + "_schema").schema
        

    def getInitOutputDatasets(self):
        # Customize init output dataset retrieval for PipelineTask.
        return {"outputSchema": SourceCatalog(self.schema)}

    def getSchemaCatalogs(self):
        # Customize schema dataset retrieval for CmdLineTask
        return {self.config.output: SourceCatalog(self.schema)}

    def _getConfigName(self):
        # Config writing with CmdLineTask is disabled for this class.
        return None

    def _getMetadataName(self):
        # Metadata writing with CmdLineTask is disabled for this class.
        return None

    def runDataRef(self, patchRefList):
        """Run this task via CmdLineTask and Gen2 Butler.
        Parameters
        ----------
        patchRefList : `list` of `lsst.daf.persistence.ButlerDataRef`
            A list of DataRefs for all filters in a single patch.
        """
        #import pdb;pdb.set_trace()
        images = {}
        replacers = {}
        mergedDataId = {"tract": patchRefList[0].dataId["tract"],
                        "patch": patchRefList[0].dataId["patch"]}
        butler = patchRefList[0].butlerSubset.butler
        ref = butler.get(self.config.ref, dataId=mergedDataId)
        imageId = butler.get("deepMergedCoaddId", dataId=mergedDataId)
        for patchRef in patchRefList:
            filt = getShortFilterName(patchRef.dataId["filter"])
            images[filt] = patchRef.get(self.config.images)

            fpCat = patchRef.get(self.config.deblendCatalog)
            footprints = {rec.getId(): (rec.getParent(), rec.getFootprint()) for rec in fpCat}
            replacers[filt] = NoiseReplacer(self.config.deblendReplacer, exposure=images[filt],
                                            footprints=footprints, exposureId=imageId)
        results = self.run(images, ref, imageId=imageId, replacers=replacers)
        butler.put(results.output, self.config.output, dataId=mergedDataId)

    def defineSchema(self, refSchema):
        """Return the Schema for the output catalog.
        This may add or modify self.
        Parameters
        ----------
        refSchema : `lsst.afw.table.Schema`
            Schema of the input reference catalog.
        Returns
        -------
        outputSchema : `lsst.afw.table.Schema`
            Schema of the output catalog.
        """
        raise NotImplementedError("Must be implemented by derived classes.")

    def run(self, images, ref):
        """Process coadds from all bands for a single patch.
        This method should not add or modify self.
        Parameters
        ----------
        images : `dict` of `lsst.afw.image.ExposureF`
            Coadd images and associated metadata, keyed by filter name.
        ref : `lsst.afw.table.SourceCatalog`
            A catalog with one record for each object, containing "best"
            measurements across all bands.
        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Struct with (at least) an `output` attribute that is a catalog
            to be written as ``self.config.output``.
        """
        raise NotImplementedError("Must be implemented by derived classes.")
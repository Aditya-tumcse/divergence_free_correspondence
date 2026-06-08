from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMakeDeps


class DivFreeCorrespondenceConanfile(ConanFile):
    settings = "os", "compiler", "build_type", "arch"

    def requirements(self):
        self.requires("pcl/1.13.1")
        self.requires("ceres-solver/2.2.0")
        self.requires("eigen/3.4.0", override=True)

    def build_requirements(self):
        self.test_requires("catch2/3.5.2")

    def configure(self):
        # Disable heavy optional PCL deps not needed for this project
        self.options["pcl"].with_opengl = False
        self.options["pcl"].with_png = False
        self.options["pcl"].with_pcap = False
        self.options["ceres-solver"].use_cuda = False

    def generate(self):
        CMakeDeps(self).generate()
        CMakeToolchain(self).generate()

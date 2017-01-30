import abc

class Designer:
	__metaclass__ = abc.ABCMeta

	@abc.abstractmethod
	def design(self, spec):
		'''Base design abstraction'''
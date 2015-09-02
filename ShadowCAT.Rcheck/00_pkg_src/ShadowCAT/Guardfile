module ::Guard
  class RPackage < Plugin
    # Initializes a Guard plugin.
    # Don't do any work here, especially as Guard plugins get initialized even
    # if they are not in an active group!
    #
    # @param [Hash] options the custom Guard plugin options
    # @option options [Array<Guard::Watcher>] watchers the Guard plugin file watchers
    # @option options [Symbol] group the group this Guard plugin belongs to
    # @option options [Boolean] any_return allow any object to be returned from a watcher
    #
    def initialize(options = {})
      super
    end

    # Called once when Guard starts. Please override initialize method to init stuff.
    #
    # @raise [:task_has_failed] when start has failed
    # @return [Object] the task result
    #
    def start
      ::Guard::UI.info('Starting R-Package guard')
    end

    # Called when just `enter` is pressed
    # This method should be principally used for long action like running all specs/tests/...
    #
    # @raise [:task_has_failed] when run_all has failed
    # @return [Object] the task result
    #
    def run_all
      test_package
      build_and_check
    end

    # Default behaviour on file(s) changes that the Guard plugin watches.
    # @param [Array<String>] paths the changes files or paths
    # @raise [:task_has_failed] when run_on_change has failed
    # @return [Object] the task result
    #
    def run_on_changes(_paths)
      test_package
    end

    # Called on file(s) additions that the Guard plugin watches.
    #
    # @param [Array<String>] paths the changes files or paths
    # @raise [:task_has_failed] when run_on_additions has failed
    # @return [Object] the task result
    #
    def run_on_additions(_paths)
      test_package
    end

    # Called on file(s) modifications that the Guard plugin watches.
    #
    # @param [Array<String>] paths the changes files or paths
    # @raise [:task_has_failed] when run_on_modifications has failed
    # @return [Object] the task result
    #
    def run_on_modifications(_paths)
      test_package
    end

    # Called on file(s) removals that the Guard plugin watches.
    #
    # @param [Array<String>] paths the changes files or paths
    # @raise [:task_has_failed] when run_on_removals has failed
    # @return [Object] the task result
    #
    def run_on_removals(_paths)
      test_package
    end

    private

    def build_and_check
      build
      check
    end

    def build
      ::Guard::UI.info('Building R-Package')
      system('R CMD build .')
    end

    def check
      ::Guard::UI.info('Checking R-Package')
      system('R CMD check ShadowCAT_0.1.tar.gz')
    end

    def test_package
      ::Guard::UI.info('Running R-Package tests')
      system("R --no-save --quiet -e 'devtools::test()'")
    end
  end
end


# clearing :on

guard 'r-package' do

  watch(%r{R/(.+)\.R$})
  watch(%r{tests/testthat/(.*)\.R})
end
